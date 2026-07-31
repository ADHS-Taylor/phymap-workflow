library(dplyr)
library(tidyr)

# Mock location lookup
location_lookup <- data.frame(
  location = c("New York", "Boston"),
  lat = c(40.7, 42.3),
  long = c(-74.0, -71.1)
)

# Return 0-row data frame to simulate empty lookup / error / NA
snap_to_nearest_empty <- function(q_lat, q_lon, lookup) {
  lookup[integer(0), ]
}

df <- data.frame(
  startLat = c(40.8),
  startLon = c(-74.1),
  endLat = c(42.2),
  endLon = c(-71.0),
  startYear = c(2025.1)
)

df_with_info <- df %>%
  rowwise() %>%
  mutate(
    start_info = list(snap_to_nearest_empty(startLat, startLon, location_lookup)),
    end_info   = list(snap_to_nearest_empty(endLat, endLon, location_lookup))
  ) %>%
  ungroup()

print("Before unnest with 0-row df:")
print(df_with_info)
print(str(df_with_info$start_info))

tryCatch({
  df_unnested <- df_with_info %>%
    unnest_wider(start_info, names_sep = "_") %>%
    unnest_wider(end_info, names_sep = "_")
  print("After unnest:")
  print(df_unnested)
  print(colnames(df_unnested))
}, error = function(e) {
  print(paste("Error occurred:", e$message))
})
