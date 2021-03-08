#' consumption_from_xml: calculate consumption from pairs of plates
#'
#'
#' @param imp_xml: An internal xml document imported via xml2::read_xml().
#'
#' vol: The initial volume dispensed into the wells of the microplate.
#' Default vol = 10.
#'
#' matrix_type: A character variable indicating the type of built-in
#' layout used to organize the 96-well culture plate. Default matrix_type = "custom".
#'
#' For custom matrix layouts, matrix_type = "custom".
#' layout_matrix A matrix containing a user-provided layour matrix for
#' sample in the 96-well culture plate. Default layout_matrix = NULL.
#'
#' Example layout matrix with C = CSB, E = Evaporation, M = Male, F = Female
#' des_matrix = matrix(c("C", "E", "M", "M", "M", "M", "M", "F", "F", "F", "F", "F",
#'                       "C", "E", "M", "M", "M", "M", "M", "F", "F", "F", "F", "F",
#'                       "C", "E", "M", "M", "M", "M", "M", "F", "F", "F", "F", "F",
#'                       "C", "E", "M", "M", "M", "M", "M", "F", "F", "F", "F", "F",
#'                       "C", "E", "M", "M", "M", "M", "M", "F", "F", "F", "F", "F",
#'                       "C", "E", "M", "M", "M", "M", "M", "F", "F", "F", "F", "F",
#'                       "C", "E", "M", "M", "M", "M", "M", "F", "F", "F", "F", "F",
#'                       "C", "E", "M", "M", "M", "M", "M", "F", "F", "F", "F", "F"),
#'                     nrow = 8,
#'                     ncol = 12,
#'                     byrow = TRUE)
#'
#' calc_duration: A logical variable indicating whether to calculate the time
#' interval between spectrophotometer readings. Default calc_duration = FALSE.
#'
#' soln_index: A 4 element list indicating the location of solutions
#' in the 4 wells from which the flies drink.
#' Default soln_index = c("S", "N", "N", "C"), with S = Sucrose, N = Nothing,
#' and C = Cocaine. The order indicates the following layout:
#'  S N
#'  N C
#'
#' @keywords read, absorbance, layout, consumption
#'
#' @export
#'
#' @examples
#' imp_xml = xml2::read_xml(xml_doc_path)
#'
#' layout_matrix = matrix(c("C", "E", "M", "M", "M", "M", "M", "F", "F", "F", "F", "F",
#'                          "C", "E", "M", "M", "M", "M", "M", "F", "F", "F", "F", "F",
#'                          "C", "E", "M", "M", "M", "M", "M", "F", "F", "F", "F", "F",
#'                          "C", "E", "M", "M", "M", "M", "M", "F", "F", "F", "F", "F",
#'                          "C", "E", "M", "M", "M", "M", "M", "F", "F", "F", "F", "F",
#'                          "C", "E", "M", "M", "M", "M", "M", "F", "F", "F", "F", "F",
#'                          "C", "E", "M", "M", "M", "M", "M", "F", "F", "F", "F", "F",
#'                          "C", "E", "M", "M", "M", "M", "M", "F", "F", "F", "F", "F"),
#'                          nrow = 8,
#'                          ncol = 12,
#'                          byrow = TRUE)
#'
#'
#' con_data = consumption_from_xml(imp_xml)
#'
#' Output:
#' ID               Proj    Line  Sex   Meth  Group Rep   Solution      Cons(uL)
#' DGRP_0112_2MFA   DGRP    0112   M    2MFA    M    1        S        1.28602090
#' DGRP_0112_2MFA   DGRP    0112   M    2MFA    M    1        C        0.84100040
#' DGRP_0112_2MFA   DGRP    0112   M    2MFA    M    2        S        0.60428090
#' DGRP_0112_2MFA   DGRP    0112   M    2MFA    M    2        C        1.02102090
#' DGRP_0112_2MFA   DGRP    0112   M    2MFA    F    1        S        1.81433180
#' DGRP_0112_2MFA   DGRP    0112   M    2MFA    F    1        C        2.17056277
#' DGRP_0112_2MFA   DGRP    0112   M    2MFA    F    2        S        1.18248175
#' DGRP_0112_2MFA   DGRP    0112   M    2MFA    F    2        C        1.38156028
#'
consumption_from_xml = function(imp_xml,
                                vol = 10,
                                matrix_type = "custom",
                                layout_matrix = NULL,
                                calc_duration = FALSE,
                                soln_index = c("S", "N", "N", "C"),
                                sample_ID) {

  plate_array = read_plates(imp_xml,
                            calc_duration = calc_duration,
                            soln_index = soln_index,
                            layout_matrix = layout_matrix,
                            sample_ID = sample_ID)

  consumption = consumption_formula(vol, plate_array[,,1], plate_array[,,2])

  parse_IDs(str_replace(dimnames(plate_array[,,])[[3]][1], " \\d", ""))

  consumption = cbind(parse_IDs(str_replace(dimnames(plate_array[,,])[[3]][1], " \\d", "")),
                      consumption)

  colnames(consumption) = c("ID","Block", "Line", "Analyst", "Group", "Well.ID", "Row", "Col", "Solution", "Cons")

  return(consumption)
}



#' read_plates: creates a 3 dimensional array containing the plate reads
#' for each sample in a [plate row, plate column, plate ID] format
#'
#'
#' @param imp_xml: An internal xml document imported via xml2::read_xml().
#'
#'
#' @keywords read, absorbance, threshold
#'
#' @export
#'
#' @examples
#' imp_xml = xml2::read_xml(xml_doc_path)
#'
#' abs_reads = read_plates(imp_xml)
#'
#' print(abs_reads)
#'
#' Output:
#'      [,1]   [,2]   [,3]   [,4]   ...[,12]
#' [1,] 0.2774 0.1596 0.1390 0.4679
#' [2,] 0.1667 0.1481 0.3128 0.1542
#' [3,] 0.1526 0.2950 2.2933 0.1468
#' [4,] 0.1509 0.1556 0.1655 0.2997
#' ...
#' [8,]
read_plates = function(imp_xml, calc_duration, soln_index, layout_matrix, sample_ID) {
  nodes = xml_path(xml_find_all(imp_xml, str_glue("/Experiment/PlateSections/PlateSection[@Name=\"{sample_ID}\"]")))

    # Iterate over the nodesets for each ID indicated by i
    for (j in 1:length(nodes)){

      # Retrieve and reorganize data from the node
      abs = list_data(imp_xml = imp_xml,plate_node = nodes[j], soln_index = soln_index, layout_matrix = layout_matrix)

      # Initialize or add to the array holding the plates and the indices
      if (j == 1){
        sample_plates = array(abs, dim = c(nrow(abs), ncol(abs), 1))
        mat_index = paste(sample_ID, " ", j)
      } else {
        sample_plates = abind(sample_plates,abs, along=3)
        mat_index = rbind(mat_index,paste(sample_ID, " ", j))
      }
    }

  # Set the indices for the matrices in the sample plates array
  dimnames(sample_plates)[[3]] = mat_index
  return(sample_plates)
}



#' sample_list: generates a simplified sample list containing unique IDs
#' and optional duration between reads.
#'
#'
#' @param imp_xml: An internal xml document imported via xml2::read_xml().
#'
#' calc_duration: A logical variable indicating whether to return
#' estimated consumption duration. Default calc_duration = FALSE.
#'
#' @keywords consumption time, elapsed time
#'
#' @export
#'
#' @examples
#' imp_xml = xml2::read_xml(xml_doc_path)
#'
#' # Get all reads for double checking
#' all_reads = import_plates(imp_xml)
#'
#' print(all_reads)
#'
#' Output:
#' ID               Project   Line    Assay   ReadTime
#' HARD_0000_SPEC    HARD     0000    SPEC    2020-07-02 11:54:00
#' HARD_0000_SPEC    HARD     0000    SPEC    2020-07-03 08:49:00
#' HARD_TAPE_SPEC    HARD     TAPE    SPEC    2020-07-02 14:10:00
#' HARD_TAPE_SPEC    HARD     TAPE    SPEC    2020-07-03 10:53:00
#'
#' # Calculate consumption times from matching IDs
#' consumption_times = sample_list(imp_xml, calc_duration = TRUE)
#'
#' print(consumption_times)
#'
#' ID               Project  Line       Assay   Duration(hr)
#' HARD_0000_SPEC    HARD     0000      SPEC    20.9166666666667
#' HARD_TAPE_SPEC    HARD     TAPE      SPEC    20.7166666666667
#'
sample_list = function(imp_xml, calc_duration) {
  # Import data from internal xml document
  all_reads = import_plates(imp_xml)

  # Identify rows with matching sample IDs
  # Create dataframe to store time calculations for each ID
  summary.data = data.frame(ID = unique(all_reads$ID))

  if (calc_duration == TRUE) {
    comp_df = con_time(all_reads, summary.data)
  } else {
    comp_df = parse_IDs(summary.data$ID)
  }

  return(comp_df)
}



#' import_plates: generates a full readout of all plate reads in an
#' xml document. Calls parse_IDs and parse_time functions to reformat
#' the xml attributes.
#'
#'
#' @param imp_xml An internal xml document imported via xml2::read_xml().
#'
#' @keywords import, plates, read information
#'
#' @export
#'
#' @examples
#' imp_xml = xml2::read_xml(xml_doc_path)
#' all_reads = import_plates(imp_xml)
#' print(all_reads)
#'
#' Output:
#' ID               Project  Line       Assay    ReadTime
#' HARD_0420_BLZR    HARD    0420       BLZR     2020-07-02 11:54:00
#' TIPS_8008_SPEC    TIPS    8008       SPEC     2020-07-03 08:49:00
#' TORK_TAPE_SPEC    TORK    TAPE       SPEC     2020-07-02 14:10:00
#' HARD_BAPE_SPEC    HARD    BAPE       SPEC     2020-07-03 10:53:00
import_plates = function(imp_xml) {
  # Identify all read plates
  plates = xml_children(imp_xml)

  # Retrieve read information from individual plates
  read.info = as.data.frame(xml_attrs(xml_children(plates)))

  # Format read information for easier handling
  read.info = apply(read.info,2, unlist)
  colnames(read.info) = NULL
  read.info = as.data.frame(t(read.info))

  # Transpose SampleIDs and append to read.info dataframe
  read.info = cbind(parse_IDs(read.info$Name), select(read.info, -c("Name", "InstrumentInfo")))

  # Parse read time for subsequent calculations
  read.info$ReadTime = parse_time(read.info$ReadTime)

  return(as.data.frame(read.info))
}



#' parse_IDs: converts a vector of sample IDS into a matrix of
#' factors for each respective ID.
#'
#' @param id: A vector containing character strings of sample IDs
#' with each factor separated by a character. sep The character
#' separating factors within the sample ID. Default sep = "_".
#'
#' @keywords IDs, samples, factors
#'
#' @export
#'
#' @examples
#' sampleID = c("DGRP_0120_MUSH",
#'              "DGRP_5430_BLIP")
#' factors = parse_IDs(sampleID)
#' print(factors)
#'
#' Output:
#'      [,1]                [,2]    [,3]    [,4]
#' [1,] "DGRP_0120_M_MUSH", "DGRP", "0120", "MUSH"
#  [2,] "DGRP_5430_F_BLIP", "DGRP", "5430", "BLIP"
parse_IDs = function(id, sep = "_") {
  id.list = data.frame(NULL)

  # Iterates over IDs in id.list to split and combine the IDs and
  # individual factors
  for (i in seq(1,length(id))) {
    temp = unlist(strsplit(id[i],split = sep))
    temp = append(temp, id[i], after=0)
    id.list = rbind(id.list, temp)
  }

  colnames(id.list) = c("ID", "Project", "Line", "Assay")
  return(id.list)
}




#' parse_time: converts date-time to lubridate format for easier
#' time calculations
#'
#' @param date_time_string A character string containing the
#' date and time as exported from Molecular Devices iD5 plate reader.
#' i.e. HH:MM am/pm M/D/YYYY
#'
#' @keywords time, date
#'
#' @export
#'
#' @examples
#' ReadTime = "2:54 PM 7/2/2020"
#' time = parse_time(ReadTime)
#' print(time)
#'
#' Output: 2020-07-02 14:54:00
parse_time = function(date_time_string) {
  time = readr::parse_datetime(date_time_string,"%I%.%M %p %m/%d/%Y")
  time = lubridate::force_tz(time, tzone = "US/Eastern")
  return(time)
}



#' list_data: extract readings from xml document and returns a tidy-formatted list.
#'
#'
#' @param imp_xml: An internal xml document imported via xml2::read_xml().
#'
#' plate_node: The location of the node containing the plate reads that
#' are being re-formatted. Note: For clarification concerning navigating
#' xml documents, search for XPATH tutorials.
#'
#' soln_index: A 4 element list indicating the location of solutions
#' in the 4 wells from which the flies drink.
#' Default soln_index = c("S", "N", "N", "C") which corresponds to:
#'  S N
#'  N C
#'
#' @keywords
#'
#' @export
#'
#' @examples
#' imp_xml = xml2::read_xml(xml_doc_path)
#' plate_node = "/Experiment/PlateSections[2]/PlateSection"
#'
#' # Get a 32x48 matrix of absorbance readings
#' abs_1536 = make_1536(imp_xml, plate_node)
#'
#' print(abs_1536)
#'
#' Output:
#'      [,1]   [,2]   [,3]   [,4]   ...[,48]
#' [1,] 0.2774 0.1596 0.1390 0.4679
#' [2,] 0.1667 0.1481 0.3128 0.1542
#' [3,] 0.1526 0.2950 2.2933 0.1468
#' [4,] 0.1509 0.1556 0.1655 0.2997
#' ...
#' [32,]
list_data = function(imp_xml,plate_node, soln_index = c("S", "N", "N", "C"), layout_matrix) {
  # Iterate over the plate by rows
  for (i in 1:8) {
    # Iterate over each row by columns
    for (j in 1:12) {
      # Read values from the "RawData" node
      abs = xml_text(xml_find_first(imp_xml,str_glue(plate_node,"/Wavelengths/Wavelength/Wells/Well[{j+(12*(i-1))}]/RawData")))
      loc = xml_attr(xml_find_first(imp_xml,str_glue(plate_node,"/Wavelengths/Wavelength/Wells/Well[{j+(12*(i-1))}]")), "Name")
      row = xml_attr(xml_find_first(imp_xml,str_glue(plate_node,"/Wavelengths/Wavelength/Wells/Well[{j+(12*(i-1))}]")), "Row")
      col = xml_attr(xml_find_first(imp_xml,str_glue(plate_node,"/Wavelengths/Wavelength/Wells/Well[{j+(12*(i-1))}]")), "Col")

      # Process data from node (split string, convert to numeric, unlist)
      abs = strsplit(abs, split = " ")
      abs = unlist(lapply(abs, as.numeric))

      if (i == 1 & j == 1) {
        temp_frame = cbind(layout_matrix[i,j], loc, row, col, soln_index, abs)
        colnames(temp_frame) = c("Group", "Well.ID", "Row", "Col", "Solution", "Abs")
      } else {
        temp_frame = rbind(temp_frame,cbind(layout_matrix[i,j], loc, row, col, soln_index, abs))
      }
    }

  }
  return(temp_frame)
}



#' consumption_formula: formula used to calculate consumption
#'
#'
#' @param vol The initial volume dispensed into the wells of the microplate.
#' Default vol = 10.
#'
#' plate1 The pre-exposure plate in the calculation.
#'
#' plate2 The post-exposure plate in the calculation.
#'
#' @keywords formula, calculate
#'
#' @export
consumption_formula = function(vol, plate1, plate2) {
  # plate1 - pre-exposure plate reading
  # plate2 - post-exposure plate reading
  # Positive = decreased absorbance
  # Negative = increased absorbance
  plate1 = as.data.frame(cbind(plate1, plate2[,"Abs"]))

  colnames(plate1) = c("Group", "Well.ID", "Row", "Col", "Solution", "Abs1", "Abs2")

  plate1$Abs1 = as.numeric(plate1$Abs1)
  plate1$Abs2 = as.numeric(plate1$Abs2)

  plate1 = mutate(plate1, Cons = vol*((Abs1-Abs2)/Abs1))

  plate1 = select(plate1, -c("Abs1", "Abs2"))

  return(plate1)
}





#' get_stats: apply a standard set of statistical tests to the data
#'
#' @param con_list A list of consumption values belonging to a single
#' group for analysis.
#' remove_outliers A logical variables specifying whether to remove
#' outliers from the data via IQR.
#'
#' @keywords statistics, normality, descriptive
#'
#' @export
get_stats = function(x, remove_outliers = FALSE) {
  if (remove_outliers == FALSE) {
    norm_test = shapiro.test(x)
    ttest = t.test(x, mu = 0, alternative = "two.sided")
    q = quantile(x)
    upper = q[4]+1.5*IQR(x)
    lower = q[2]-1.5*IQR(x)
    outlier_range = paste0("(",round(lower,3), ")-(", round(upper,3), ")")
    outlier_loc = x > upper | x < lower
    outlier_count = length(x[outlier_loc])
    inlier_count = length(x[!outlier_loc])
    stats = rbind(median(x), mean(x), var(x), sd(x), sd(x)/sqrt(length(x)), sd(x)/mean(x), norm_test$p.value, ttest$p.value)
    rownames(stats) = c("Median", "Mean", "Var", "SD", "SEM", "Coef. of Variation", "Normality (p>0.05)", "t test vs 0")
    colnames(stats) = deparse(substitute(x))

  } else if (remove_outliers == TRUE) {
    q = quantile(x)
    upper = q[4]+1.5*IQR(x)
    lower = q[2]-1.5*IQR(x)
    outlier_loc = x > upper | x < lower
    outlier_count = length(x[outlier_loc])
    inlier_count = length(x[!outlier_loc])

    y = x[!outlier_loc]

    norm_test = shapiro.test(y)
    ttest = t.test(y, mu = 0, alternative = "two.sided")
    stats = rbind(median(y), mean(y), var(y), sd(y), sd(y)/sqrt(length(y)), sd(y)/mean(y), norm_test$p.value, ttest$p.value, upper, lower, paste0(outlier_count,"/",inlier_count))
    rownames(stats) = c("Median", "Mean", "Var", "SD", "SEM", "Coef. of Variation", "Normality (p>0.05)", "t test vs 0", "Upper Cutoff (+1.5IQR)", "Lower Cutoff (-1.5IQR)", "Outliers (Out/In)")
    colnames(stats) = deparse(substitute(x))
  }
  return(stats)
}



#' con_time: estimates the consumption time.
#'
#'
#' @param
#'
#' @keywords consumption time, elapsed time
#'
#' @export
#'
#' @examples
#' imp_xml = xml2::read_xml(xml_doc_path)
#'
#' # Get all reads for double checking
#' all_reads = import_plates(imp_xml)
#'
#' print(all_reads)
#'
#' Output:
#' ID                  Project  Line   Sex  Assay   ReadTime
#' HARD_0000_0_SPEC    HARD     0000   0    SPEC    2020-07-02 11:54:00
#' HARD_0000_0_SPEC    HARD     0000   0    SPEC    2020-07-03 08:49:00
#' HARD_TAPE_0_SPEC    HARD     TAPE   0    SPEC    2020-07-02 14:10:00
#' HARD_TAPE_0_SPEC    HARD     TAPE   0    SPEC    2020-07-03 10:53:00
#'
#' # Calculate consumption times from matching IDs
#' consumption_times = con_time(imp_xml)
#'
#' print(consumption_times)
#'
#' ID                  Project  Line   Sex  Assay   Duration(hr)
#' HARD_0000_0_SPEC    HARD     0000   0    SPEC    20.9166666666667
#' HARD_TAPE_0_SPEC    HARD     TAPE   0    SPEC    20.7166666666667
#'
con_time = function(all_reads, summary.data) {
  # Expand the logical vectors to make them pairwise comparisons
  exp_matr = get_logical_comps(all_reads, summary.data)

  # Calculate the time interval for matching IDs in logical matrix
  for (j in seq(8,ncol(exp_matr))){
    pairs = exp_matr[,j]

    temp = exp_matr[pairs,]

    temp = arrange(temp, desc(ReadTime))

    temp_time = lubridate::interval(start = temp[2,"ReadTime"], end = temp[1,"ReadTime"])

    temp_duration = lubridate::int_length(temp_time)/3600

    temp_vec = data.frame(exp_matr[min(which(pairs == TRUE)),1:6],
                          temp_duration)

    temp_vec$ReadNo = paste(exp_matr[pairs,"ReadNo"],collapse = " vs ")

    colnames(temp_vec) = c("ID", "Project", "Line", "Sex", "Assay", "Comparison", "Duration(hr)")

    if (j == 8) {
      comp_df = temp_vec
    } else {
      comp_df = rbind(comp_df, temp_vec)
    }
  }


  return(comp_df)
}



#' get_logical_comps:
#'
#'
#' @param
#'
#' @keywords
#'
#' @export
#'
#' @examples
get_logical_comps = function(all_reads, summary.data) {
  # Generate logical matrix showing which rows share the same ID
  indices = sapply(all_reads$ID,function(x) {x == summary.data[,"ID"]})
  colnames(indices) = NULL
  indices = t(indices)

  read_ind = replicate(nrow(indices),0)

  read_ind = as.vector(apply(indices, 2, function(x) { read_ind[x] = seq(1, sum(x))}))

  # Split logical vector into independent vectors each with only
  # a single value
  # Iterate over each full logical vector
  for (i in seq(1,ncol(indices))){
    vec = indices[,i]
    init = TRUE

    while (sum(vec) > 0){
      new_vec = logical(length(vec))

      loc = min(which(vec == TRUE))
      vec[loc] = FALSE
      new_vec[loc] = TRUE

      if (init){
        temp_matr = as.matrix(new_vec)
        init = FALSE
      } else {
        temp_matr = cbind(temp_matr, as.matrix(new_vec))
      }
    }

    if (ncol(temp_matr) < 2){
      temp_pairs = temp_matr
    } else{
      pairs = combn(ncol(temp_matr),2)
      temp_pairs = apply(pairs, 2, function(x) {rowSums(temp_matr[,x])})
    }


    if (i == 1) {
      exp_matr = as.matrix(temp_pairs)
    } else {
      exp_matr = cbind(exp_matr, temp_pairs)
    }

  }

  exp_matr = sapply(as.data.frame(exp_matr), as.logical)

  all_reads = add_column(all_reads,ReadNo = read_ind, .after = "Assay")

  exp_matr = cbind(all_reads, exp_matr)

  return(exp_matr)
}



#' calc_preference: calculate the preference in 2 choice MFA using
#' dataframe returned from comsumption_from_xml()
#'
#'
#' @param data_list: dataframe generated using consumption_from_xml()
#'
#' @keywords
#'
#' @export
#'
#' @examples
#'
calc_preference = function(data_list, equation = "difference") {
  grps = unique(data_list$Group)
  for (i in 1:length(grps)){
    # Generate logical array showing which rows share the group
    grp_idx = sapply(data_list$Group, function(x) {x == grps[i]})

    # Generate a sample list that's a subset of the group
    sample_list = unlist(unique(data_list[grp_idx,"ID"]))

    # Iterate over each sample ID
    for (j in 1:length(sample_list)){
      # Generate a logical array showing which rows share sample IDs
      sample_idx = sapply(data_list$ID, function(x) {x == sample_list[j]})

      # Generate a list containing the number of replicates in the group & sample
      rep_list = unlist(unique(data_list[grp_idx & sample_idx,"Well.ID"]))

      # Iterate over the number of replicates
      for (k in 1:length(rep_list)){
        # Generate a logical array showing which rows have the same replicate number
        rep_idx = sapply(data_list$Well.ID, function(x) {x == rep_list[k]})

        # Isolate the two rows containing the same group label, ID, and replicate number
        rep_data = data_list[grp_idx & sample_idx & rep_idx, ]

        if (nrow(rep_data) != 2) {
          error_name = data_list[grp_idx & sample_idx & rep_idx, c("ID", "Group", "Well.ID", "Solution")]
          warning(str_glue("Error in data calculation: incorrect number of readings for: /n/n {error_name}"))
          next
        }

        # Isolate the consumption values for sucrose and cocaine
        sucrose = rep_data[str_detect(rep_data$Solution,"S|Sucrose"), "Cons"]
        cocaine = rep_data[str_detect(rep_data$Solution,"C|Cocaine"), "Cons"]

        # Calculate cocaine preference
        if (equation == "ratio") {
          pref = cocaine/sucrose
        } else if (equation == "difference") {
          pref = (cocaine - sucrose)/(cocaine + sucrose)
        } else if( equation == "sum") {
          pref = cocaine + sucrose
        }


        # Create a new dataframe to hold the preference value
        temp_df = select(rep_data[1,],-c("Solution", "Cons"))
        Pref = data.frame(Pref = pref)
        colnames(Pref) = "Pref"
        temp_df = cbind(temp_df, Pref)

        # Initialize or build the new dataframe
        if (i == 1 & j ==1 & k == 1) {
          pref_df = temp_df
        } else {
          pref_df = rbind(pref_df, temp_df)
        }
      }
    }
  }
  return(pref_df)
}



#' evap_adjust:
#'
#'
#' @param data_list: dataframe generated using consumption_from_xml()
#'
#' evap_label: character used to indicate evaporation wells in the layout
#' matrix. Default evap_labe = "E".
#'
#' @keywords
#'
#' @export
#'
#' @examples
#'
evap_adjust = function(data_list, evap_label = "E"){
  data_list = as.data.frame(data_list)
  # Isolate out unique sample IDs (i.e. plates)
  sample_list = unique(data_list$ID)

  # Iterate over sample IDs
  for (i in 1:length(sample_list)){
    # Find rows containing the current sample ID
    sample_idx = sapply(data_list$ID, function(x) {x == sample_list[i]})

    # Isolate the rows for the current sample ID
    temp_df = data_list[sample_idx,]

    # Find rows containing evaporation wells indicated by evap_label
    evap_idx = str_detect(temp_df$Group, pattern = evap_label)

    # Subtract median evaporation value from all values
    for (j in unique(temp_df$Solution)) {
      soln_idx = str_detect(temp_df$Solution, j)
      temp_df$Cons = temp_df$Cons - median(temp_df[evap_idx & soln_idx,"Cons"])
    }

    # Drop the evaporation rows
    temp_df = temp_df[!evap_idx,]

    # Initialize or grow dataframe containing new data
    if (i == 1) {
      adj_data = temp_df
    } else {
      adj_data = rbind(adj_data, temp_df)
    }
  }

  adj_data = as_tibble(adj_data)

  return(adj_data)
}





#' csb_adjust:
#'
#'
#' @param data_list: dataframe generated using consumption_from_xml()
#'
#' evap_label: character used to indicate evaporation wells in the layout
#' matrix. Default evap_label = "E".
#'
#' @keywords
#'
#' @export
#'
#' @examples
#'
csb_adjust = function(data_list, csb_label = "C"){
  # Isolate out unique sample IDs (i.e. plates)
  sample_list = unique(data_list$ID)

  population_median_csb = median(data_list[data_list$Group == csb_label,]$Cons)

  # Iterate over sample IDs
  for (i in seq(1,length(sample_list))){
    # Find rows containing the current sample ID
    sample_idx = sapply(data_list$ID, function(x) {x == sample_list[i]})

    # Isolate the rows for the current sample ID
    temp_df = data_list[sample_idx,]

    # Find rows containing CSB wells indicated by csb_label
    csb_idx = grep(csb_label, temp_df$Group)

    # Add difference between between population CSB median and plate CSB median to the plate
    temp_df$Cons = temp_df$Cons+ (population_median_csb - median(temp_df[csb_idx,"Cons"]))

    # Drop the evaporation rows
    temp_df = temp_df[-csb_idx,]

    # Initialize or grow dataframe containing new data
    if (i == 1) {
      adj_data = temp_df
    } else {
      adj_data = rbind(adj_data, temp_df)
    }
  }

  return(adj_data)
}



#' link_xml_to_layout:
#'
#'
#' @param data_list: dataframe generated using consumption_from_xml()
#'
#' evap_label: character used to indicate evaporation wells in the layout
#' matrix. Default evap_label = "E".
#'
#' @keywords
#'
#' @export
#'
#' @examples
#'
link_xml_to_layout = function(layout_dir, xml_dir) {
  # generate list of sample IDs contained in plate layouts (master list)
  layout_list = data.frame(layout_path = file.path(layout_dir,list.files(layout_dir, pattern = ".csv$")))
  layout_list$XML_path = NA

  # Generate list of XML files
  xml_paths = file.path(xml_dir,list.files(xml_dir, pattern = ".xml$"))

  for (i in xml_paths) {
    temp_xml = xml2::read_xml(i)
    xml_IDs = xml_attr(xml_children(xml_children(temp_xml)), "Name")
    link_idx = str_replace(basename(layout_list$layout_path), ".csv", "") %in% xml_IDs

    layout_list[link_idx, "XML_path"] = i

  }
  rm(temp_xml)

  na_idx = is.na(layout_list$XML_path)

  dropped_paths = basename(layout_list[na_idx, "layout_path"])

  if (sum(na_idx)) warning(str_glue("Unlinked plate layouts in directory: {dropped_paths} \n\n"))

  return(layout_list)
}




















