import sys
import time
import pandas as pd
import pyfaidx

def main():
    time_stamp = time.time()
    script_name = sys.argv[0]
    ref_genome = sys.argv[1]
    input_file = sys.argv[2]
    output_file = sys.argv[3] + "_" + str(time_stamp)

    print("Getting guides for base editing using {}".format(script_name))
    print("for the variants from file {}".format(input_file))
    print("and writing to files {}.csv and {}_filtered.csv.".format(output_file, output_file))
    print("Using this path to your reference genome: {}".format(ref_genome))


    ref_genome_pyfaidx = pyfaidx.Fasta(ref_genome)
    variant_list = pd.read_csv(input_file, sep=',')["variant"]

    all_variant = []
    all_editable = []
    all_target_seq_ref = []
    all_target_seq = []
    all_target_base_ref = []
    all_target_base = []
    all_possible_guides_with_pam = []
    all_edit_strings = []
    all_edit_pos_strings = []
    all_possible_guides = []
    all_possible_pams = []
    all_rev_com = []

    for variant in variant_list:
        chrom, position, alt = variant.split("_")
        chrom = chrom.replace("chr", "")
        position = int(position) - 1

        length = 20

        position_edit_window_start = 4
        position_edit_window_end = 8

        PAM = "NGG" # or "NG"
        PAM_length = len(PAM) # 2 or 3

        total_length = length + PAM_length # 22 or 23
        length_edit_window = position_edit_window_end - (position_edit_window_start - 1)

        bases_before_ew = position_edit_window_start - 1
        bases_after_ew = length - position_edit_window_end
        bases_after_ew_with_PAM = total_length - position_edit_window_end

        bases_before_variant = bases_before_ew + length_edit_window - 1 # - 1 since variant always need to stay in the edit window
        bases_after_variant_with_PAM = bases_after_ew_with_PAM + length_edit_window - 1 # - 1 since variant always need to stay in the edit window

        target_base_ref = str(ref_genome_pyfaidx[chrom][position])

        if (target_base_ref == "A" and alt == "G") or (target_base_ref == "C" and alt == "T"):
            editable = True
            rev_com = False
            target_base = target_base_ref

            target_seq_ref = str(ref_genome_pyfaidx[chrom][position-bases_before_variant:position+bases_after_variant_with_PAM+1]) # + 1 since last position is excluded
            target_seq = target_seq_ref

        elif (target_base_ref == "G" and alt == "A") or (target_base_ref == "T" and alt == "C"):
            editable = True
            rev_com = True
            target_base = target_base_ref.replace('A', '*').replace('T', 'A').replace('*', 'T').replace('C', '&').replace('G', 'C').replace('&', 'G')[::-1]

            target_seq_ref = str(ref_genome_pyfaidx[chrom][position-bases_after_variant_with_PAM:position+bases_before_variant+1]) # + 1 since last position is excluded
            target_seq = target_seq_ref.replace('A', '*').replace('T', 'A').replace('*', 'T').replace('C', '&').replace('G', 'C').replace('&', 'G')[::-1]

        else:
            editable = False


        if editable:
            possible_guides_with_pam = []
            possible_guides = []
            possible_pams = []
            edit_strings = []
            edit_pos_strings = []

            count_of_N_in_PAM = PAM.count("N")
            count_of_G_in_PAM = PAM_length - count_of_N_in_PAM

            start_of_PAM_search = bases_after_variant_with_PAM - bases_after_ew - count_of_N_in_PAM # the N of NGG/NG is irrelevant for the search here
            end_of_PAM_search = count_of_G_in_PAM - 1 # the N of NGG/NG is irrelevant for the search here; stop at beginning of last PAM
            length_of_PAM_search = count_of_G_in_PAM # the N of NGG/NG is irrelevant for the search here

            start_before_hit = length + count_of_N_in_PAM
            end_after_hit_with_PAM = count_of_G_in_PAM
            end_after_hit_guide_only = count_of_N_in_PAM

            variant_position = bases_before_variant

            for i in range(len(target_seq)-start_of_PAM_search,
                        len(target_seq)-end_of_PAM_search):
                if target_seq[i: i+length_of_PAM_search] == PAM.replace("N", ""):
                    possible_guide_with_pam = target_seq[i-start_before_hit:i+end_after_hit_with_PAM]
                    possible_guides_with_pam.append(possible_guide_with_pam)
                    possible_guide = target_seq[i-start_before_hit:i-end_after_hit_guide_only]
                    possible_guides.append(possible_guide)
                    possible_pam = target_seq[i-count_of_N_in_PAM:i+count_of_G_in_PAM]
                    possible_pams.append(possible_pam)
                    edit_string = "-" * length
                    edit_string = list(edit_string)
                    edit_pos_string = "-" * length
                    edit_pos_string = list(edit_pos_string)
                    for j in range(position_edit_window_start - 1, position_edit_window_end):
                        if possible_guide[j] == target_base and j == variant_position:
                            edit_string[j] = "V"
                            edit_pos_string[j] = str(j+1)
                        elif possible_guide[j] == target_base:
                            edit_string[j] = "*"
                            edit_pos_string[j] = str(j+1)
                    edit_string = "".join(edit_string)
                    edit_strings.append(edit_string)
                    edit_pos_string = "".join(edit_pos_string)
                    edit_pos_strings.append(edit_pos_string)
                variant_position -= 1

            all_variant.append(variant)
            all_editable.append(editable)
            all_target_seq_ref.append(target_seq_ref)
            all_target_seq.append(target_seq)
            all_target_base_ref.append(target_base_ref)
            all_target_base.append(target_base)
            all_possible_guides_with_pam.append(possible_guides_with_pam)
            all_edit_strings.append(edit_strings)
            all_edit_pos_strings.append(edit_pos_strings)
            all_possible_guides.append(possible_guides)
            all_possible_pams.append(possible_pams)
            all_rev_com.append("{}".format(rev_com))

        else:
            all_variant.append(variant)
            all_editable.append(editable)
            all_target_seq_ref.append("")
            all_target_seq.append("")
            all_target_base_ref.append("")
            all_target_base.append("")
            all_possible_guides_with_pam.append("")
            all_edit_strings.append("")
            all_edit_pos_strings.append("")
            all_possible_guides.append("")
            all_possible_pams.append("")
            all_rev_com.append("{}".format(""))

    sgrnas = pd.DataFrame({"all_variant": all_variant,
                        "all_editable": all_editable,
                        "all_target_seq_ref": all_target_seq_ref,
                        "all_target_seq": all_target_seq,
                        "all_target_base_ref": all_target_base_ref,
                        "all_target_base": all_target_base,
                        "all_possible_guides_with_pam": all_possible_guides_with_pam,
                        "all_edit_strings": all_edit_strings,
                        "all_edit_pos_strings": all_edit_pos_strings,
                        "all_possible_guides": all_possible_guides,
                        "all_possible_pams": all_possible_pams,
                        "all_rev_com": all_rev_com})

    sgrnas.to_csv(output_file + ".csv", sep=';', index=False)

    sgrnas_filtered = sgrnas[sgrnas.astype(str)["all_possible_guides"] != "[]"]

    sgrnas_filtered.to_csv(output_file + "_filtered.csv", sep=';', index=False)

if __name__ == "__main__":
    main()
