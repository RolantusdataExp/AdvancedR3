#' Descriptive statistics for lipidomics data
#'
#' @param data
#'
#' @return A data.frame/tibble
#'
#'
#'
descriptive_stats <- function(data) {
    data %>%
        dplyr::group_by(metabolite) %>%
        dplyr::summarise(dplyr::across(value, list(mean = mean, sd = sd))) %>%
        dplyr::mutate(dplyr::across(tidyselect::where(is.numeric), ~ round(.x, digits = 1)))
}

#' Plot histogram distribution
#'
#' @param Lipidomics object
#'
#' @return Plot by ggplot2
#'
#'
#'
plot_distributions = function(data) {
    data %>%
        ggplot2::ggplot( ggplot2::aes(x = value)) + ggplot2::geom_histogram() +            ggplot2::facet_wrap(ggplot2::vars(metabolite), scales = "free")
}

#' Change to snakecase a given column!
#'
#' @param dataframe
#' @param cols - coloumn in dataframe
#'
#' @return
#'
#'
#' @examples
column_values_to_snake_case <- function(data, cols) {
    data %>%
        dplyr::mutate(dplyr::across({{ cols }}, snakecase::to_snake_case))
}

#' Function to create a pivote-wider dataformat
#'
#' @param data
#'
#' @return pivot-wider
#' @export
#'
#' @examples
metabolites_to_wider = function(data) {
    data %>%
        pivot_wider(names_from = metabolite, values_from = value, values_fn = mean, names_prefix = "metabolite_")
}

#' A transformation recipe to pre-process the data.
#'
#' @param data The lipidomics dataset.
#' @param metabolite_variable The column of the metabolite variable.
#'
#' @return
#'
create_recipe_spec <- function(data, metabolite_variable) {
    recipes::recipe(data) %>%
        recipes::update_role({{ metabolite_variable }}, age, gender, new_role = "predictor") %>%
        recipes::update_role(class, new_role = "outcome") %>%
        recipes::step_normalize(tidyselect::starts_with("metabolite_"))
}

#' Create model
#'
#' @param model_specs is the model specs
#' @param recipe_specs recipe specs
#'
#' @return workflow object
#' @export
#'
#' @examples
create_model_workflow = function(model_specs, recipe_specs) {
    workflows::workflow() %>%
        workflows::add_model(model_specs) %>%
        workflows::add_recipe(recipe_specs)
}

#' Tidy model output
#'
#' @param workflow_fitted_model
#'
#' @return Something super exciting!
#' @export
#'
#' @examples
tidy_model_output = function(workflow_fitted_model) {
    workflow_fitted_model %>%
        workflows::extract_fit_parsnip() %>%
        broom::tidy(exponentiate = TRUE)
}

#' Convert the long form dataset into a list of wide form
#'
#' @param Data lipidomics
#'
#' @return List of dataframes
#'
#'
#'
split_by_metabolite <- function(data) {
    data %>%
        column_values_to_snake_case(metabolite) %>%
        dplyr::group_split(metabolite) %>%
        purrr::map(metabolites_to_wider)
}

#' Generate the results of a model
#'
#' @param data The lipidomics dataset.
#'
#' @return A data frame.
#'
generate_model_results <- function(data) {
    create_model_workflow(
        parsnip::logistic_reg() %>%
            parsnip::set_engine("glm"),
        data %>%
            create_recipe_spec(tidyselect::starts_with("metabolite_"))
    ) %>%
        parsnip::fit(data) %>%
        tidy_model_output()
}

#' Title calculates the model estimates as well as the code that adds the original metabolite names into functions
#'
#' @param model_results glm model results
#' @param data lipidomics dataset
#'
#' @return
#' @export a dataframe
#'
#' @examples
add_original_metabolite_names = function(model_results, data) {
    data %>%
        dplyr::select(metabolite) %>%
        dplyr::mutate(term = metabolite) %>%
        column_values_to_snake_case(term) %>%
        dplyr::mutate(term = stringr::str_c("metabolite_", term)) %>%
        dplyr::distinct(term, metabolite) %>%
        right_join(model_results, by = "term")
}


#' Title Calculate the estimates for the model for each metabolite.
#'
#' @param data lipidomics
#'
#' @return
#' @export
#'
#' @examples
calculate_estimates <- function(data){
    data %>%
        split_by_metabolite() %>%
        purrr::map(generate_model_results) %>%
        purrr::list_rbind() %>%
        dplyr::filter(stringr::str_detect(term, "metabolite_")) %>%
        add_original_metabolite_names(data)
}

#' Title Function for plotting estimates of all metabolites
#'
#' @param results estimates of model
#'
#' @return object
#'
#'
#'
plot_estimates <- function(results) {
    results %>%
        ggplot2::ggplot(aes(
            x = estimate,
            y = metabolite,
            xmin = estimate - std.error,
            xmax = estimate + std.error
        )) + ggplot2::geom_pointrange() + ggplot2::coord_fixed(xlim = c(0,5))
}

