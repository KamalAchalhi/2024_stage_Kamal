"""
Récuperation des fichiers nécessaire (Atlas, images sujet anatomiques), recalage rigide et calcul de similarité
Pour trouver le meilleur template/atlas pour chaque sujet. Et retourner la listes de bons atlas et la transformations nécessaires.
"""

import os
import time
import numpy as np
import tools as tls
import matplotlib.pyplot as plt
def etape1(nom_general_sujet, path_pattern, path_ouput, Nom_Atlas_anatomic_Pattern, path_des_atlas, file_transfo_direc, file_transfo_inv):
    """
    :param nom_general_sujet: Pattern du nom d'image d'un sujet anatomique données
    :param path_pattern: Pattern du chemin pour parvenir à chaque sujet
    :param path_ouput: Chemin pour sauver chaque images après traitement
    :param Nom_Atlas_anatomic_Pattern: Nom d'un atlas anatomique quelquonque
    :param path_des_atlas: Repertoire contenant tout les atlas anatomiques
    :param file_transfo_direc: Repertoire pour sauvegarder les transformations de recalage (sujet to atlas)
    :param file_transfo_inv: Repertoire pour sauvegarder les transformations de recalage (atlas to sujet)
    :return:
            list_atlas_meilleur : listes path du meilleurs atlas pour chaque sujet
            tab_img_sujet_rot : Listes path pour sujet après traitmeent swapping
            list_tranf_direc : Listes Path pour chaque transformation directe (sujet i to his best atlas)
            list_tranf_inv :  Listes Path pour chaque transformation inverse (best atlas to his sujet i)
    """

    ###COMMENTER LE MAXIMUM DE LIGNE POUR RENDRE CLAIRE
    debut = time.time()
    files_atlas = tls.Parcours_dossier_only_data_match(path_des_atlas, Nom_Atlas_anatomic_Pattern)   #
    tab_path_sujet = tls.recup_les_sujets(nom_general_sujet, pattern_sous_repertoire_by_sujet=path_pattern)  #
    print(tab_path_sujet)
    list_path_sujet_rot = [tls.creation_PATH_pour_fichier_swaper(sujet_path, path_ouput) for sujet_path in tab_path_sujet] #
    for path_sujet_rot, path_sujet in zip(list_path_sujet_rot, tab_path_sujet):
        tls.SWAP_COPY_INFO_SAVE(path_sujet, path_sujet_rot)
    tab_repertoire, tab_img_sujet = tls.path_abs_sujet_to_fichier_repertorie_sujet(list_path_sujet_rot) #
    print(tab_img_sujet)

    criteres = ['MattesMutualInformation']
    list_atlas_meilleur = list()
    list_tranf_direc = []
    list_tranf_inv = []
    for sujet, repertoire in zip(tab_img_sujet, tab_repertoire):
        bon_atlas, path_trf_direct, path_trf_inv = tls.recup_bon_atlas_avc_transfos(files_atlas, criteres, path_des_atlas, sujet, repertoire, "Rigid", "linear", file_transfo_direc, file_transfo_inv)
        list_atlas_meilleur.append(bon_atlas)
        list_tranf_direc.append(path_trf_direct)
        list_tranf_inv.append(path_trf_inv)
        print(f"l'atlas qui maximise l'information mutuel est : {bon_atlas} pour {sujet}\n")
    fin = time.time()
    tps_excecution = fin - debut
    print(f"le temps d'exécution du programme est : {tps_excecution} secondes")
    return list_atlas_meilleur, list_path_sujet_rot, list_tranf_direc, list_tranf_inv


if __name__ == "__main__":
    file_transfo_direc = r'/envau/work/meca/users/2024_Kamal/output/output_script1'
    file_transfo_inv = r'/envau/work/meca/users/2024_Kamal/output/output_script2'
    path_variables = "/envau/work/meca/users/2024_Kamal/2024_stage_Kamal/variables"
    nom_general_sujet = r'^sub-00\d+\_ses-00\d+\_acq-haste_rec-nesvor_desc-aligned_T2w.nii.gz'
    path_pattern = r'/envau/work/meca/users/2024_Kamal/real_data/lastest_nesvor/sub-00\d+\/ses-00\d+\/haste/default_reconst'
    all_sujets_path = "/envau/work/meca/users/2024_Kamal/real_data/lastest_nesvor"
    path_output = "/envau/work/meca/users/2024_Kamal/output/output_script1"
    Nom_Atlas_anatomic_Pattern = r'^STA\d+\.nii.gz'
    path_des_atlas = "/envau/work/meca/users/2024_Kamal/Sym_Hemi_atlas/Fetal_atlas_gholipour/T2"
    path_output_avec_mask = "/envau/work/meca/users/2024_Kamal/output/Output_script1_avec_mask"
    # chemin_vers_mask = np.load(os.path.join(path_variables,"chemin_vers_mask.npy"))
    list_atlas_meilleur, list_path_sujet_rot, list_tranf_direc, list_tranf_inv = etape1(nom_general_sujet, path_pattern, path_output, Nom_Atlas_anatomic_Pattern, path_des_atlas, file_transfo_direc, file_transfo_inv)

    print(list_atlas_meilleur, list_tranf_inv)
    np.save(os.path.join(path_variables, "list_atlas_meilleur.npy"), list_atlas_meilleur, allow_pickle='False')
    np.save(os.path.join(path_variables, "list_path_sujet_rot.npy"), list_path_sujet_rot, allow_pickle='False')
    np.save(os.path.join(path_variables, "list_tranf_direc.npy"), list_tranf_direc, allow_pickle='False')
    np.save(os.path.join(path_variables, "list_tranf_inv.npy"), list_tranf_inv, allow_pickle='False')
