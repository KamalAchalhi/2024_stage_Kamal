import os
import re
import time
import numpy as np
import pandas as pd
import nibabel as nib
import ants
import matplotlib.pyplot as plt
# import matplotlib
# matplotlib.use('Qt5Agg')
import tools as tls



def etape3(SUB_rec_by_Atlas_PATH, path_output):
    """

    :param SUB_rec_by_Atlas_PATH:
    :param path_output:
    :return:
    """
    debut = time.time()
    list_repertoire, list_image_sub_recal = tls.path_abs_sujet_to_fichier_repertorie_sujet(SUB_rec_by_Atlas_PATH)
    list_path_threshold = []
    for image_sub_recal, repertoire in zip(list_image_sub_recal, list_repertoire):
        Image_recal = nib.load(os.path.join(repertoire, image_sub_recal))
        Image_recal_array = Image_recal.get_fdata()
        Image_recal_array_out = np.zeros_like(Image_recal_array)
        Image_recal_array_out[Image_recal_array == 2] = 10
        #Image_recal_array[Image_recal_array != 10] = 0
        path_image_threshold = tls.creation_chemin_nom_img(path_output, image_sub_recal, "seg_L_only_x10.nii.gz")
        list_path_threshold.append(path_image_threshold)
        Image_recal_threshold = nib.Nifti1Image(Image_recal_array_out, Image_recal.affine, Image_recal.header) # nous devons copier les infos géometrique pour aligner des img qui sont deja ds meme espace
        nib.save(Image_recal_threshold, path_image_threshold)
    fin = time.time()
    tps_excecution = fin - debut
    print(f"le temps d'exécution du programme est : {tps_excecution} secondes")
    return list_path_threshold


if __name__ == "__main__":
    path_variables = "/envau/work/meca/users/2024_Kamal/2024_stage_Kamal/variables"
    path_output = "/envau/work/meca/users/2024_Kamal/output/output_script3"
    SUB_rec_by_Atlas_PATH = np.load(os.path.join(path_variables,"SUB_rec_by_Atlas_PATH.npy"))
    list_path_threshold = etape3(SUB_rec_by_Atlas_PATH, path_output)
    np.save(os.path.join(path_variables, "list_path_threshold.npy"), list_path_threshold, allow_pickle='False')



