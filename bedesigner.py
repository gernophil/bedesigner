import os
import sys
import argparse
import time
import pandas as pd
import pyfaidx

def main():

    parser = argparse.ArgumentParser(description='bedesigner',
                                     prog='python bedesigner.py')

    required = parser.add_argument_group("required arguments")
    required.add_argument('-i', '--input', help='Path to the input file',
                          metavar='<PATH>', required=True)
    required.add_argument('-g', '--ref-genome', help='Path to the reference genome',
                          metavar='<PATH>', required=True)
    required.add_argument('-o', '--output', help='Path to the output file',
                          metavar='<PATH>', required=True)

    parser.add_argument('-p', '--pam-site', help='Sequence of the PAM site',
                        choices=['NGG', 'NG'], default='NGG', type=str)
    parser.add_argument('-s', '--window-start', help='Starting position of editing window',
                        default=4, type=int)
    parser.add_argument('-e', '--window-end', help='End position of editing window',
                        default=8, type=int)
    parser.add_argument('-l', '--guide-length', help='Length of guide',
                        default=20, type=int)
    parser.add_argument('-n', '--ignore-string', help='Substring to be ignored in the variant string',
                        type=str)
    parser.add_argument('-b', '--base-editor', help='Base editor to design guides for',
                        choices=['both', 'ABE', 'CBE'], default='both', type=str)
    parser.add_argument('-a', '--all-possible', help='Design for all ALTs',
                        action='store_true')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()

    args = parser.parse_args()

    if not os.path.isfile(args.input):
        sys.exit('Input file not found!')

    if not os.path.isfile(args.ref_genome):
        sys.exit('Reference genome not found!')


    time_stamp = time.time()
    script_name = sys.argv[0]
    refgenome = args.ref_genome
    input_file = args.input
    pamsite = args.pam_site
    edit_window_start = args.window_start
    edit_window_end = args.window_end
    guidelength = args.guide_length
    ignorestring = args.ignore_string
    output_file = args.output + "_p" + pamsite + "_e" + \
                  str(edit_window_start) + "-" + str(edit_window_end) + \
                  "_l" + str(guidelength) # + "_" + str(time_stamp)

    baseeditor = args.base_editor
    if baseeditor == "both":
        ABE = True
        CBE = True
    elif baseeditor == "ABE":
        ABE = True
        CBE = False
    elif baseeditor == "CBE":
        ABE = False
        CBE = True

    print("Getting guides for base editing using {}".format(script_name))
    print("for the variants from file {}".format(input_file))
    print("and writing to files {}.csv and {}_filtered.csv.".format(output_file, output_file))
    print("Using this path to your reference genome: {}".format(refgenome))


    ref_genome_pyfaidx = pyfaidx.Fasta(refgenome)
    variant_list = pd.read_csv(input_file, sep=',')["variant"]

    all_variant = []
    all_editable = []
    all_be_strings = []
    all_original_alt = []
    all_target_seq_ref = []
    all_target_seq_ref_match = []
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
        variant_coords = variant.split("_")
        if len(variant_coords) == 3:
            chrom, position, alt = variant_coords
            alt = alt.upper()
        elif len(variant_coords) == 4:
            chrom, position, ref, alt = variant_coords
            alt = alt.upper()
            ref = ref.upper()
        else:
            raise ValueError("Your variants are improperly formatted. Please use [CHROM]_[POS]_[ALT] or [CHROM]_[POS]_[REF]_[ALT].")
        if ignorestring:
            chrom = chrom.replace(ignorestring, "")
        position = int(position) - 1

        length = guidelength

        position_edit_window_start = edit_window_start
        position_edit_window_end = edit_window_end

        PAM = pamsite
        PAM_length = len(PAM) # 2 or 3

        total_length = length + PAM_length
        length_edit_window = position_edit_window_end - (position_edit_window_start - 1)

        bases_before_ew = position_edit_window_start - 1
        bases_after_ew = length - position_edit_window_end
        bases_after_ew_with_PAM = total_length - position_edit_window_end

        bases_before_variant = bases_before_ew + length_edit_window - 1 # - 1 since variant always need to stay in the edit window
        bases_after_variant_with_PAM = bases_after_ew_with_PAM + length_edit_window - 1 # - 1 since variant always need to stay in the edit window

        target_base_ref = str(ref_genome_pyfaidx[chrom][position])
        if ref:
            if ref == target_base_ref:
                target_seq_ref_match = "ref_match"
            else:
                target_seq_ref_match = "REF_MISMATCH"
        else:
            target_seq_ref_match = "no_ref_input"

        if (ABE and target_base_ref == "A" and alt == "G") or \
           (ABE and target_base_ref == "A" and args.all_possible) or \
           (CBE and target_base_ref == "C" and alt == "T") or \
           (CBE and target_base_ref == "C" and args.all_possible):
            editable = True
            rev_com = False
            if target_base_ref == "A":
                be_string = "ABE"
                if alt == "G":
                    original_alt = True
                else:
                    original_alt = False
            elif target_base_ref == "C":
                be_string = "CBE"
                if alt == "T":
                    original_alt = True
                else:
                    original_alt = False
            target_base = target_base_ref

            target_seq_ref = str(ref_genome_pyfaidx[chrom][position-bases_before_variant:position+bases_after_variant_with_PAM+1]) # + 1 since last position is excluded
            target_seq = target_seq_ref

        elif (ABE and target_base_ref == "T" and alt == "C") or \
             (ABE and target_base_ref == "T" and args.all_possible) or \
             (CBE and target_base_ref == "G" and alt == "A") or \
             (CBE and target_base_ref == "G" and args.all_possible):
            editable = True
            rev_com = True
            if target_base_ref == "T":
                be_string = "ABE"
                if alt == "C":
                    original_alt = True
                else:
                    original_alt = False
            elif target_base_ref == "G":
                be_string = "CBE"
                if alt == "A":
                    original_alt = True
                else:
                    original_alt = False
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
            all_be_strings.append(be_string)
            all_original_alt.append(original_alt)
            all_target_seq_ref.append(target_seq_ref)
            all_target_seq_ref_match.append(target_seq_ref_match)
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
            all_be_strings.append("")
            all_original_alt.append("")
            all_target_seq_ref.append("")
            all_target_seq_ref_match.append(target_seq_ref_match)
            all_target_seq.append("")
            all_target_base_ref.append("")
            all_target_base.append("")
            all_possible_guides_with_pam.append("")
            all_edit_strings.append("")
            all_edit_pos_strings.append("")
            all_possible_guides.append("")
            all_possible_pams.append("")
            all_rev_com.append("{}".format(""))

    sgrnas = pd.DataFrame({"variant": all_variant,
                           "editable": all_editable,
                           "base_editor": all_be_strings,
                           "originally_intended_ALT": all_original_alt,
                           "target_seq_reference_genome": all_target_seq_ref,
                           "ref_match": all_target_seq_ref_match,
                           "targeted_seq_by_guide": all_target_seq,
                           "target_base_reference_genome": all_target_base_ref,
                           "targeted_base_by_guide": all_target_base,
                           "all_possible_guides_with_pam": all_possible_guides_with_pam,
                           "all_off-target_bases_by_guide": all_edit_strings,
                           "all_edited_positions_by_guide": all_edit_pos_strings,
                           "all_possible_guides": all_possible_guides,
                           "all_possible_pams": all_possible_pams,
                           "guide_target_opposite_strand": all_rev_com})

    sgrnas.to_csv(output_file + ".csv", sep=';', index=False)

    sgrnas_filtered = sgrnas[(sgrnas.astype(str)["all_possible_guides"] != "[]") &
                             (sgrnas.astype(str)["all_possible_guides"] != "")]

    sgrnas_filtered.to_csv(output_file + "_filtered.csv", sep=';', index=False)

if __name__ == "__main__":
    main()
