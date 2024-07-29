"""
Application de la transformations de recalage inverse sur un atlas de segmentation d'hemisphère LR vers le sujet associès (dont il est le plus similaire)
"""
import os
import time
import numpy as np
import nibabel as nib
import ants
import tools as tls


def etape2(path_des_atlas_hemi_seg, list_atlas_meilleur, list_path_sujet_rot, list_tranf_inv, path_output_repertoire):
    """

    :param path_des_atlas_hemi_seg: Chemins des atlas segmenté par hemisphères RL
    :param list_atlas_meilleur: Listes des chemins pour chaque atlas mieux adaptés à un sujet données
    :param list_path_sujet_rot: Listes chemin vers repertoire où on était sauvé les iamges de sujet anatomiques après "swaping" du script1
    :param list_tranf_inv: Listes des chemins vers les transformations de recalage inverse obtenue entre chaque sujet et son meilleur atlas
    :param path_output_repertoire: Chemin de sauvegarde des images après traitement réalisé par script2
    :return: AtlasLR_rec_dans_sub_space : Atlas segmenté LR recalé dans l'espace sujet associès
    """
    debut = time.time()

    tab_repertoire, tab_img_sujet = tls.path_abs_sujet_to_fichier_repertorie_sujet(list_path_sujet_rot)  #Separation path repertoire du nom de l'image
    tls.creation_data_frame_sujet_by_best_atlas(tab_img_sujet, list_atlas_meilleur)  #Creation tableau panda contenant pour chaque sujet son meilleur atlas associés
    les_atlas_binary = [] #liste qui recuperera les path des atlas hemisphère binaire
    list_num = tls.extraction_numero_atlas(list_atlas_meilleur) # récupération pour chaque atlas adapté du numero/age de l'atlas dans une liste
    for num in list_num:  #parcours de chaque age
        les_atlas_binary.append(f'STA{num}_all_reg_LR_dilM.nii.gz') # rassemble format nom des atlas hemisphere binaire et de l'age recupere, obtention liste des atlas binaire gauche droite adapté
    AtlasRL_rec_dans_sub_space= []  # liste pour conserver les chemins des atlas segmentés LR recaler dans l'espace du sujet
    for sujet, repertoire,  atlas_binar, warp in zip(tab_img_sujet, tab_repertoire, les_atlas_binary, list_tranf_inv): # parcours les sujet swaper, les atlas binar RL et les transfo inverse associées
        Sujet_fixe = ants.image_read(os.path.join(repertoire, sujet)) #recup l'image sujet swap en format ants
        Atlas_binary = ants.image_read(os.path.join(path_des_atlas_hemi_seg, atlas_binar))  # recup atlas binar RL grâce au nom formé ligne29 et path athlas binar RL
        transfo = warp + '_Inverse_0GenericAffine.mat' # étant donné que ma fct enregistre les fonction de recalage dans la list sans prendre en compte ce sufixe, il doit être rajouté
        Atlas_binary_warped = ants.apply_transforms(Sujet_fixe, Atlas_binary,  transformlist=transfo, direction = 'fwd', interpolator="nearestNeighbor") #Application de la trf inverse sur l'atlas binar RL adapté vers le ième sujet
        path_atlas_binary_warped = tls.creation_chemin_nom_img(path_output_repertoire, sujet, atlas_binar)  # creation de path de l'atlas binar RL apres son recalage dans espace sujet
        AtlasRL_rec_dans_sub_space.append(path_atlas_binary_warped)  # rajoute à liste des path des atlas RL apres rec
        ants.image_write(Atlas_binary_warped, path_atlas_binary_warped)  # Enregistre l'image obtenue dans le path crée
    fin = time.time()
    tps_excecution = fin - debut
    print(f"le temps d'exécution du programme est : {tps_excecution} secondes")
    return AtlasRL_rec_dans_sub_space


if __name__ == "__main__":
    nom_general_sujet = r'^sub-00\d+\_ses-00\d+\_acq-haste_rec-nesvor_desc-aligned_T2w.nii.gz'  #f ormat nom sujet anatomique
    path_pattern = r'/envau/work/meca/users/2024_Kamal/real_data/lastest_nesvor/sub-00\d+\/ses-00\d+\/haste/default_reconst'  # format path repertoire sujet anatomique
    path_des_atlas_hemi_seg = r'/envau/work/meca/users/2024_Kamal/Sym_Hemi_atlas'  # path des atlas binaire par hemisphère
    path_variables = "/envau/work/meca/users/2024_Kamal/2024_stage_Kamal/variables"
    path_repertoire_sujet_rot = "/envau/work/meca/users/2024_Kamal/output/output_script1"
    path_output_repertoire = "/envau/work/meca/users/2024_Kamal/output/output_script2"  # path du repertoire pour sauver les atlas RL apres recalage
    list_atlas_meilleur = np.load(os.path.join(path_variables, "list_atlas_meilleur.npy"))   # recup list des meilleurs atlas
    list_path_sujet_rot = np.load(os.path.join(path_variables, "list_path_sujet_rot.npy"))
    list_tranf_direc =  np.load(os.path.join(path_variables, "list_tranf_direc.npy"))
    list_tranf_inv = np.load(os.path.join(path_variables, "list_tranf_inv.npy"))

    AtlasRL_rec_dans_sub_space = etape2(path_des_atlas_hemi_seg, list_atlas_meilleur, list_path_sujet_rot, list_tranf_inv, path_output_repertoire)
    np.save(os.path.join(path_variables, "AtlasRL_rec_dans_sub_space.npy"), AtlasRL_rec_dans_sub_space, allow_pickle='False')