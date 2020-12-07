## Utility functions used during model training and prediction

#' Record model predictions and associated data
#' 
#' @param trained_model a trained model (of any class for which a predict() method exists)
#' @param new_data data frame for which predictions should be generated
#' @param data_name a name for this dataset (e.g. "test data")
#' @param data_cols vector of column names, specifying columns from `new_data` which should
#'                  be recorded alongside predictions
#' @param iteration unique identifier for the current iteration
#' 
record_predictions <- function(trained_model, new_data, data_name, data_cols, iteration) {
  new_data %>% 
    mutate(iteration = iteration,
           dataset = data_name,
           prediction = predict(trained_model, newdata = ., na.action = na.pass),
           prob = predict(trained_model, newdata = ., type = 'prob', na.action = na.pass)[["True"]]) %>% 
    select(.data$iteration, .data$dataset, all_of(data_cols), .data$prediction, .data$prob)
}
