args <- commandArgs(trailingOnly = TRUE)

library(tidyverse)

bedesigner_filt <- read_delim(args[1])

bedesigner_long <-
  bedesigner_filt |>
  separate_longer_delim(c(all_possible_guides_with_pam,
                          `all_off-target_bases_by_guide`,
                          all_edited_positions_by_guide,
                          all_possible_guides,
                          all_possible_pams), ",") |>
  mutate(all_possible_guides_with_pam = gsub("[", "",
                                             gsub("'", "",
                                                  gsub("]", "",
                                                       gsub(" ", "",
                                                            all_possible_guides_with_pam,
                                                            fixed = TRUE),
                                                       fixed = TRUE),
                                                  fixed = TRUE),
                                             fixed = TRUE)) |>
  mutate(`all_off-target_bases_by_guide` = gsub("[", "",
                                                gsub("'", "",
                                                     gsub("]", "",
                                                          gsub(" ", "",
                                                               `all_off-target_bases_by_guide`,
                                                               fixed = TRUE),
                                                          fixed = TRUE),
                                                     fixed = TRUE),
                                                fixed = TRUE)) |>
  mutate(all_edited_positions_by_guide = gsub("[", "",
                                              gsub("'", "",
                                                   gsub("]", "",
                                                        gsub(" ", "",
                                                             all_edited_positions_by_guide,
                                                             fixed = TRUE),
                                                        fixed = TRUE),
                                                   fixed = TRUE),
                                              fixed = TRUE)) |>
  mutate(all_possible_guides = gsub("[", "",
                                    gsub("'", "",
                                         gsub("]", "",
                                              gsub(" ", "",
                                                   all_possible_guides,
                                                   fixed = TRUE),
                                              fixed = TRUE),
                                         fixed = TRUE),
                                    fixed = TRUE)) |>
  mutate(all_possible_pams = gsub("[", "",
                                  gsub("'", "",
                                       gsub("]", "",
                                            gsub(" ", "",
                                                 all_possible_pams,
                                                 fixed = TRUE),
                                            fixed = TRUE),
                                       fixed = TRUE),
                                  fixed = TRUE)) |>
  group_by(variant) |>
  mutate(guide_name = paste0(variant, "_", row_number())) |>
  ungroup()

guides_only <-
  bedesigner_long |>
  select(guide_name,
         all_possible_guides) |>
  rename(guide_sequence = all_possible_guides)

write_tsv(bedesigner_long, args[2])
write_tsv(guides_only, args[3])
