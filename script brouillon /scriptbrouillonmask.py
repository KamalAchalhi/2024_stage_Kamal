import os
import re
import time
import numpy as np
import pandas as pd
import nibabel as nib
import ants
import matplotlib.pyplot as plt
import tools as tls

def etape1bis(all_sujets_path, path_output, nom_general_sujet ):
    debut = time.time()
    list_path_img_segmente = tls.recup_les_sujets(nom_general_sujet, repertoire_sujet_segm = all_sujets_path)
    list_path_img_segmente_rot = [tls.creation_PATH_pour_fichier_swaper(sujet_path, path_output) for sujet_path in list_path_img_segmente]
    for path_sujet_rot, path_sujet in zip(list_path_img_segmente_rot, list_path_img_segmente):
        tls.SWAP_COPY_INFO_SAVE(path_sujet, path_sujet_rot)
    tab_repertoire_rot, tab_img_sujet_rot = tls.path_abs_sujet_to_fichier_repertorie_sujet(list_path_img_segmente_rot)
    threshold = 1
    chemin_vers_mask = [tls.creation_chemin_nom_img(path_output, img_name, "masque_par_threshold.nii.gz")for img_name in tab_img_sujet_rot]
    for path_sujet_segm_rot, chemin_mask in zip(list_path_img_segmente_rot, chemin_vers_mask):
        Img_segm = nib.load(path_sujet_segm_rot)
        dtype_img_sujet_segm = Img_segm.get_data_dtype()
        img_sujet_segmente_array = Img_segm.get_fdata()
        img_sujet_segmente_array[img_sujet_segmente_array >= 1] = 1
        img_sujet_segmente_array[img_sujet_segmente_array == 0] = 0
        mask = img_sujet_segmente_array > 0
        # mask =mask.astype(dtype_img_sujet_segm)
        Masque_nifiti = nib.Nifti1Image(mask, Img_segm.affine, Img_segm.header)
        nib.save(Masque_nifiti,chemin_mask)
    fin = time.time()
    tps_excecution = fin - debut
    print(f"le temps d'exécution du programme est : {tps_excecution} secondes")
    return chemin_vers_mask




if __name__ == "__main__":
    nom_general_sujet = r'^sub-00\d+\_ses-00\d+\_acq-haste_rec-nesvor_desc-aligned_T2w.nii.gz'
    path_pattern =  r'/envau/work/meca/users/2024_Kamal/real_data/lastest_nesvor/sub-00\d+\/ses-00\d+\/haste/default_reconst'
    all_sujets_path = "/envau/work/meca/users/2024_Kamal/real_data/lastest_nesvor"
    path_variables = "/envau/work/meca/users/2024_Kamal/2024_stage_Kamal/variables"
    tab_img_sujet_rot = np.load(os.path.join(path_variables, "tab_img_sujet_rot.npy"))
    path_output = "/envau/work/meca/users/2024_Kamal/output/output_script1"
    path_output_avec_mask = "/envau/work/meca/users/2024_Kamal/output/Output_script1_avec_mask"
    chemin_vers_mask = etape1bis(all_sujets_path, path_output_avec_mask, nom_general_sujet)
    np.save(os.path.join(path_variables, "chemin_vers_mask.npy"), chemin_vers_mask, allow_pickle='False')