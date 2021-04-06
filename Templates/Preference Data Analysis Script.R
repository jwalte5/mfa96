# HEADER --------------------------------------------------------------
# SCRIPT NAME: Preference Data Analysis Script.R
#
# DATE CREATED: 31 March 2021
#
# AUTHOR: Joshua Walters
#
# DESCRIPTION: This script analyzes data collected using the
# Microplate Feeder Assay.
#
# NOTES
# Prior to analyzing data, prepare a working directory containing
# the following structure:
    # Working Directory/
    #                   XML/ (contains xml data files)
    #                   Plate Layouts/ (contains plate layout csv files)
    #                   Data Files/ (contains saved dataframes)
    #                   Figures and Tables/ (contains output tables and figures)


## SET DIRECTORIES AND LOAD FILES--------------------------------------

# Set working directory and locate relevant directories
work_dir = "C:/Users/jwalt/Google Drive/Work Directory/Research/Microplate Feeder Assay/ZA - Preference Coupler Preference/"
# As an alternative, the following line can be used to select
# the working directory using a window prompt
#work_dir = choose.dir() # Choose working directory

xml_dir = paste0(work_dir, "XML/")
layout_dir = paste0(work_dir, "Plate Layouts/")
data_dir = paste0(work_dir, "Data Files/")
output_dir = paste0(work_dir, "Figures and Tables/")
setwd(work_dir)

# Load MFA functions into memory by sourcing the function script
source("C:/Users/jwalt/Google Drive/Work Directory/Research/GitHub/mfa96/R/mfafunctions96.R")

# Load Libraries
library("lubridate")
library("tidyverse")
library("xml2")
library("abind")
library("ggpubr")
library("ggsci")
library("rstatix")


## IMPORT XML DATA USING MFA FUNCTIONS---------------------------------

# Create dataframe linking plate layouts to their respective XML files
data_paths = link_xml_to_layout(layout_dir = layout_dir, xml_dir = xml_dir) %>%
  drop_na()

# Iterate over each plate layout and its corresponding XML file to
# import the data using the consumption_from_xml() function
for (i in 1:nrow(data_paths)) {

  # Import the consumption dataframe from a given row of the data_paths variable
  # into a temporary dataframe
  temp_df = consumption_from_xml( imp_xml = xml2::read_xml(data_paths[i,"XML_path"]),
                                  matrix_type = "custom",
                                  layout_matrix = read.csv(data_paths[i,"layout_path"], header = FALSE, colClasses = "character"),
                                  vol = 10,
                                  calc_duration = FALSE,
                                  soln_index = c("S", "N", "N", "C"),
                                  sample_ID = basename(str_replace(data_paths[i,"layout_path"], ".csv", "")))

  # Initialize or add to the main dataframe
  if (i == 1) {
    con_df = temp_df
  } else {
    con_df = rbind(con_df, temp_df)
  }

}

# Filter out extraneous readings, such as 1536-wells with no solution (i.e. "N" indexed)
# or 96-wells marked for omission (i.e. "X")
con_df = con_df %>%
  filter(Solution != "N", Group != "X")

# Note: For readability, values in columns can be renamed using the recode() function
# e.g.:
con_df$Group = recode(con_df$Group, E = "Evaporation", F = "Female", M = "Male")
con_df$Solution = recode(con_df$Solution, S = "Sucrose", C = "Cocaine")

# For plots, variable order can be enforced using the levels attributes of the factor() function
# where the order of the levels dictates the order of the variables on the plot
# e.g.:
con_df$Group = factor(con_df$Group, levels = c("Evaporation", "Male", "Female"))
con_df$Solution = factor(con_df$Solution, levels = c("Sucrose", "Cocaine"))

## GENERATE RAW DATA SUMMARY STATISTICS AND A RAW DATA PLOT------------

# Save the dataframe as a csv
write.csv(con_df, file = file.path(data_dir, "raw_con_data.csv"))

# Generate a raw data plot and save to the output directory
jpeg(filename = paste0(output_dir, "raw_data_plot.jpg"))
ggerrorplot(con_df,
            x = "Solution",
            y = "Cons",
            fill = "Solution",
            color = "Solution",
            facet.by = "Group",
            desc_stat = "mean_se",
            size = 0.75)+
  geom_jitter(size = 1, alpha = 0.25, width = 0.2)
dev.off()

# Transform the dataframe to calculate total consumption
con_df = con_df %>%
  pivot_wider(names_from = "Solution", values_from = "Cons") %>%
  rowwise() %>%
  mutate(Total = Sucrose + Cocaine)

# Generate summary stats and save as a csv file
write.csv(data.frame(Data_Statistics = ""),file = paste0(output_dir, "stats.csv"), row.names = FALSE)
cat("\nRaw_Data_Statistics\n", file = paste0(output_dir, "stats.csv"), append = TRUE)
con_df %>%
  group_by(Group) %>%
  get_summary_stats(Sucrose, Cocaine, Total) %>%
  write.table(file = paste0(output_dir, "stats.csv"), append = TRUE, sep = ",", row.names = FALSE)


## ADJUST CONSUMPTION VALUES FOR EVAPORATION---------------------------

# Reconfigure the con_df dataframe
con_df = con_df %>%
  pivot_longer(cols = c("Sucrose", "Cocaine"), names_to = "Solution", values_to = "Cons") %>% # Change format to be compatible with evap_adjust() function
  evap_adjust(evap_label = "Evaporation") %>% # Adjust for evaporation
  pivot_wider(names_from = "Solution", values_from = "Cons") %>% # Revert dataframe back to previous form
  select(-Total) %>% # Drop previously calculated total since it used the pre-adjusted values
  rowwise() %>% # Specify to enforce element-wise column calculations
  mutate(Total = Sucrose + Cocaine) # Recalculate the total volume

## GENERATE ATTRITION STATISTICS---------------------------------------

cat("\nAttrition_Statistics", file = paste0(output_dir, "stats.csv"), append = TRUE)

drop_df = con_df %>%
  filter(Sucrose < 0 | Cocaine < 0)

con_df = con_df %>%
  filter(Sucrose >= 0 & Cocaine >= 0)

cat("\nRetained_Samples\n", file = paste0(output_dir, "stats.csv"), append = TRUE)
con_df %>%
  group_by(Group) %>%
  get_summary_stats(Sucrose, Cocaine, Total) %>%
  write.table(file = paste0(output_dir, "stats.csv"), append = TRUE, sep = ",", row.names = FALSE)

cat("\nDropped_Samples\n", file = paste0(output_dir, "stats.csv"), append = TRUE)
drop_df %>%
  group_by(Group) %>%
  get_summary_stats(Sucrose, Cocaine, Total) %>%
  write.table(file = paste0(output_dir, "stats.csv"), append = TRUE, sep = ",", row.names = FALSE)

attr_rate = nrow(drop_df)/(nrow(drop_df) + nrow(con_df))

cat(c("\n\nAttrition Rate", attr_rate), file = paste0(output_dir, "stats.csv"), sep = ",", append = TRUE)

# Save the dataframe as a csv
write.csv(con_df, file = file.path(data_dir, "adjusted_con_data.csv"))


## GENERATE EVAPORATION-ADJUSTED DATA PLOTS----------------------------

# Reconfigure dataframe to plotting-friendly format
con_df = con_df %>%
  pivot_longer(cols = c("Sucrose", "Cocaine"), values_to = "Cons", names_to = "Solution")

jpeg(filename = paste0(output_dir, "adjusted_data_plot.jpg"))
ggerrorplot(con_df,
            x = "Solution",
            y = "Cons",
            fill = "Solution",
            color = "Solution",
            facet.by = "Group",
            desc_stat = "mean_se",
            size = 0.75)+
  geom_jitter(size = 1, alpha = 0.25, width = 0.2)
dev.off()


## CALCULATE PREFERENCE AND GENERATE PREFERENCE STATISTICS AND FIGURES-

# Calculate preference
pref_df = con_df %>%
  pivot_wider(names_from = "Solution", values_from = "Cons") %>%
  rowwise() %>%
  mutate(Pref = (Cocaine - Sucrose)/Total)

pref_df = pref_df[,c("ID", "Block", "Line", "Analyst", "Group", "Well.ID", "Row", "Col", "Sucrose", "Cocaine", "Total", "Pref")]

# Save the dataframe as a csv
write.csv(pref_df, file = file.path(data_dir, "pref_data.csv"))

cat("\n\nPreference_Statistics\n", file = paste0(output_dir, "stats.csv"), append = TRUE)
pref_df %>%
  group_by(Group) %>%
  get_summary_stats(Total, Pref) %>%
  write.table(file = paste0(output_dir, "stats.csv"), append = TRUE, sep = ",", row.names = FALSE)


