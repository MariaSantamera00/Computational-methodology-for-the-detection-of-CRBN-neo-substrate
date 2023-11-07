## This script cut the different parts of the ternary complex with CRBN

load TFM_MariaSantamera/validation_set/aligned_pdb/4ci1_al.pdb; select crbn, chain B; remove not crbn; sele crbn_union, resi 321-426; remove not crbn_union; save 4ci1_crbn.pdb

reinitialize

load TFM_MariaSantamera/validation_set/aligned_pdb/5fqd_al.pdb; select crbn, chain B; select ck1a, chain C; select talido, resnam E1438; select crbn_union, crbn and resi 321-426; save 5fqd_lig_union.pdb, crbn_union or talido; save 5fqd_target.pdb, ck1a

reinitialize

load TFM_MariaSantamera/validation_set/aligned_pdb/6h0f_al.pdb; select crbn, chain B; select ikzf1, chain C; select poli, resnam Y70 and chain B; select crbn_union, crbn and resi 321-426; save 6h0f_lig_union.pdb, crbn_union or poli; save 6h0f_target.pdb, ikzf1

reinitialize

load TFM_MariaSantamera/validation_set/aligned_pdb/6h0g_al.pdb; select crbn, chain B; select zf4, chain C; select poli, resnam Y70 and chain B; select crbn_union, crbn and resi 321-426; save 6h0g_lig_union.pdb, crbn_union or poli; save 6h0g_target.pdb, zf4

reinitialize

load TFM_MariaSantamera/validation_set/aligned_pdb/5hxb_al.pdb; select crbn, chain C; select gspt1, chain A; select cc885, resnam 85C and chain C; select crbn_union, crbn and resi 321-426; save 5hxb_lig_union.pdb, crbn_union or cc885; save 5hxb_target.pdb, gspt1

reinitialize

load TFM_MariaSantamera/validation_set/aligned_pdb/7lps_al.pdb; select crbn, chain B; select ikzf2, chain C; select alv1, resnam RN9 and chain B; select crbn_union, crbn and resi 321-426; save 7lps_lig_union.pdb, crbn_union or alv1; save 7lps_target.pdb, ikzf2

reinitialize

load TFM_MariaSantamera/validation_set/aligned_pdb/6uml_al.pdb; select crbn, chain C; select sall4, chain E; select tali, resnam Y70; select crbn_union, crbn and resi 321-426; save 6uml_lig_union.pdb, crbn_union or tali; save 6uml_target.pdb, sall4

reinitialize

load TFM_MariaSantamera/validation_set/aligned_pdb/6bn7_al.pdb; select crbn, chain B; select brd4, chain C; select prota, resnam RN3; select crbn_union, crbn and resi 321-426; save 6bn7_lig_union.pdb, crbn_union or prota; save 6bn7_target.pdb, brd4

reinitialize

load TFM_MariaSantamera/validation_set/aligned_pdb/5yiz_al.pdb; select crbn, chain A; select crbn_union, crbn and resi 321-426; save 5yiz_union.pdb, crbn_union



