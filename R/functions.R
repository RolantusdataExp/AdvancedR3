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

