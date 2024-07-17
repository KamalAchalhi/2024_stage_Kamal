"""
Addition de l'atlas segmenté LR réhaussé à une image segmenté de chaque sujet associés
pour obtenir une image de sujets qui combine segmentation des zones cerebrales et des hémisphères LR.

"""
import os
import time
import numpy as np
import nibabel as nib
import tools as tls


def etape4(nom_general_sujet, all_sujets_path,path_output, tab_path_sujet, AtlasRL_rec_dans_sub_space):
    """
    :param nom_general_sujet: Pattern du nom d'image d'un sujet segmenté donnée
    :param all_sujets_path: Chemin contenant les sous-repertoires de chaque sujet segmenté
    :param path_output: Chemin de sauvegarde des images après traitement
    :param tab_img_sujet: Listes des noms pour chaque sujet anatomique(même nom que les versions segmentées)
    :param AtlasLR_rec_dans_sub_space: Listes des chemins pour les atlas segmentés par hemisphères recalés dans l'espace sujet

    """
    debut = time.time()
    print(tab_path_sujet)
    tab_repertoire, tab_img_sujet = tls.path_abs_sujet_to_fichier_repertorie_sujet(tab_path_sujet)
    print(tab_repertoire, tab_img_sujet)
    #Recuperation des images segmentés et on les swap :
    list_path_img_segmente = tls.recup_les_sujets(nom_general_sujet, repertoire_sujet_segm=all_sujets_path)
    print(list_path_img_segmente)
    list_path_img_segmente_rot = [tls.creation_PATH_pour_fichier_swaper(sujet_path, path_output) for sujet_path in list_path_img_segmente]
    for path_sujet_rot, path_sujet in zip(list_path_img_segmente_rot, list_path_img_segmente):
        tls.SWAP_COPY_INFO_SAVE(path_sujet, path_sujet_rot)
    print(AtlasRL_rec_dans_sub_space)
    
    # ON additionne apres passage en numpy array les tableau de l'image segmenté du sujet l'image de l'hemisphère droit*10
    for path_sujet_segm, AtlasRL_rec_dans_sub_space, sujet in zip(list_path_img_segmente_rot, AtlasRL_rec_dans_sub_space, tab_img_sujet):
        img_sujet_segmente = nib.load(path_sujet_segm)
        dtype_img_sujet_segm = img_sujet_segmente.get_data_dtype()
        img_sujet_segmente_array = img_sujet_segmente.get_fdata()
        AtlasLR_rec_dans_sub_space_array = nib.load(AtlasRL_rec_dans_sub_space).get_fdata()
        img_sujet_segm_binar_combined_array = img_sujet_segmente_array + 100 * AtlasLR_rec_dans_sub_space_array
        img_sujet_segm_binar_combined_array[img_sujet_segm_binar_combined_array == 100] = 0
        img_sujet_segm_binar_combined_array[img_sujet_segm_binar_combined_array == 200] = 0
        img_sujet_segm_binar_combined_array = img_sujet_segm_binar_combined_array.astype(dtype_img_sujet_segm)
        image_segm_final = nib.Nifti1Image(img_sujet_segm_binar_combined_array, img_sujet_segmente.affine, img_sujet_segmente.header)
        path_img_final = tls.creation_chemin_nom_img(path_output, sujet, "segmentation_LR.nii.gz")
        nib.save(image_segm_final, path_img_final)
    fin = time.time()
    tps_excecution = fin - debut
    print(f"le temps d'exécution du programme est : {tps_excecution} secondes")


if __name__ == "__main__":
    path_variables = "/envau/work/meca/users/2024_Kamal/2024_stage_Kamal/variables"
    nom_general_sujet = r'^sub-00\d+\_ses-00\d+\_acq-haste_rec-nesvor_desc-aligned_T2w.nii.gz'
    all_sujets_path = r"/envau/work/meca/users/2024_Kamal/real_data/lastest_nesvor/"
    path_output = "/envau/work/meca/users/2024_Kamal/output/output_script4"
    AtlasRL_rec_dans_sub_space = np.load(os.path.join(path_variables, "AtlasRL_rec_dans_sub_space.npy"))
    tab_path_sujet = np.load(os.path.join(path_variables, "list_path_sujet_rot.npy"))

    etape4(nom_general_sujet, all_sujets_path, path_output, tab_path_sujet, AtlasRL_rec_dans_sub_space )
