## Utility functions for loading model output across different folders
## - Output has the format: dataset/response_var/feature_set/...

require(readr)
require(dplyr)
require(stringr)

load_single_run <- function(dataset, response_var, run_id) {
  file.path("output", dataset, response_var, run_id, "predictions.rds") %>% 
    read_rds() %>% 
    mutate(dataset = dataset,
           response_var = str_to_sentence(response_var),
           run_id = run_id)
}

load_all_runs <- function(data_sub_set) {
  # output folder structure is: output/dataset/response/run_id/
  top_dir <- file.path("output", data_sub_set)
  response_vars <- list.dirs(top_dir, full.names = FALSE, recursive = FALSE)
  
  loaded_data <- tibble()
  
  for (resp in response_vars) {
    response_path <- file.path(top_dir, resp)
    run_ids <- list.dirs(response_path, full.names = FALSE, recursive = FALSE)
    
    new_data <- lapply(run_ids, load_single_run, 
                       dataset = data_sub_set,
                       response_var = resp) %>% 
      bind_rows()
    
    loaded_data <- bind_rows(loaded_data, new_data)
  }
  
  loaded_data
}
