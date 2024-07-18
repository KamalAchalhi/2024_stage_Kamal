import os
import re
import numpy as np
import nibabel as nib
import ants
import matplotlib.pyplot as plt
import pandas as pd
"""
je vais expliquer ce que fais chaque fonction
"""


def recup_les_sujets(nom_general_sujet, repertoire_sujet_segm = None, pattern_sous_repertoire_by_sujet = None):
    """

    :param nom_general_sujet: format nom des sujet recherché
    :param repertoire_sujet_segm: argument qui prends le path du repertoire où se trouve les sujets segemnté
    :param pattern_sous_repertoire_by_sujet: argument qui prends le path des repertoires où se trouve les sujets anatomiques
    :return:
    """
    file_paths = []
    pattern = re.compile(nom_general_sujet) #Crée un pattern correspondant au format du nom à rechercher
    if repertoire_sujet_segm :  #1er cas possible correpondant à l'une des deux utilisisation de cet algo, si l'argument feedé est le repertoire des sujet segm
        for file in os.listdir(repertoire_sujet_segm):  #On parcours chaque fichier du repertoire
             if pattern.match(file):  #s'il y a correspondance entre un fichier donnée et le pattern cherché
                    file_paths.append(os.path.join(repertoire_sujet_segm ,file))  # on rajoute à la list des path, le path complet du fichier
    elif pattern_sous_repertoire_by_sujet: #2ieme cas, si c'est l'argument du pattern des differents sous repertoires contenant les sujets qui est feedé
        path_pattern = re.compile(pattern_sous_repertoire_by_sujet)  #creation 2ieme pattern correspondant au format des paths des sujets anatomiques
        for root, _, files in os.walk("/envau/work/meca/users/2024_Kamal/real_data/lastest_nesvor/"): #On parcours chaque sous repertoire et fichiers du repertoire initial
            if path_pattern.search(root): # Si le sous repertoire correspond au format recherché alors :
                for file in files : #on y parcours tout les fichiers du repertoire correspodant
                    if pattern.match(file): #Si un fichier correspond au format de nom rechercher pour sujet anatomique, --> on l'ajoute
                        file_paths.append(os.path.join(root, file))
    return sorted(file_paths)


def copy_info_geo(path_img_input, path_img_input_copied):
    """

    :param path_img_input: chemin de l'image dont on veut copier les infos geometriques
    :param path_img_input_copied: chemin de l'image qui va recevoir les infos geometriques
    :return:
    """
    img_copied = nib.load(path_img_input_copied)  #charge l'image en format nifti
    img_input = nib.load(path_img_input).get_fdata()  #charge l'image en format numpy array necessaire pour recevoir les infos geo
    return nib.Nifti1Image(img_input, img_copied.affine, img_copied.header)  # on converti en format nifti un numpy array tout en copiant les infos geo provenant de l'image donneuse resté en format nifti


def SWAP_COPY_INFO_SAVE(path_img_input, path_img_rot_nifti):
    """

    :param path_img_input: chemin image initial avant traitement
    :param path_img_rot_nifti: chemin image apres traitement par swap pour la sauver
    :return:
    """
    img_input = nib.load(path_img_input)  # Image input à swaper en format nifti
    img_input_array = img_input.get_fdata()     # conversion en numpy array pour realiser le swap
    original_data_type = img_input.get_data_dtype()  # on recupère le data type de l'image format nifti pr pas le perdre
    img_input_array_rot = np.transpose(img_input_array, (0, 1, 2))[::1, ::1, ::-1]  # transpose np array en inversant le dernier axe
    img_input_array_rot = img_input_array_rot.astype(original_data_type)        # impose le data type initial au tableau
    img_input_rot_nifti = nib.Nifti1Image(img_input_array_rot, img_input.affine, img_input.header)  # On converti en nifti en copiant
    nib.save(img_input_rot_nifti, path_img_rot_nifti)  # On sauve l'image swapé dans le path donnée en input


def creation_PATH_pour_fichier_swaper(path_sujet, repertoire_output):
    """

    :param path_sujet: path des images initial avant traitement
    :param repertoire_output: path de sortie pour sauver les image apres traitement par swapping
    :return:
    """
    fichier = os.path.basename(path_sujet)  # on recupere uniquement le nom du fichier
    nom_initial, fin = (fichier[:-7], ".nii.gz") if fichier.endswith(".nii.gz") else os.path.splitext(fichier)  # on separe le nom en tant que tel et le suffixe (nii.gz de nifiti)
    return os.path.join(repertoire_output, f"{nom_initial}_rot{fin}")  # création d'un nouveau path vers le repertoire donné en input et
                                                                       # le nom recuperer ci dessus avec le suffixe rot pour decrire le traitement et le suffixe nifiti


def calcul_similarity_ants(img1, img2, critere, path_mask = None):
    """

    :param img1: 1er image input
    :param img2: 2ieme image input
    :param critere: metrique de similarite utitise pour calcul de similarité entre image
    :param path_mask: si un chemin pour un masque est donné en input
    :return:
    """
    return ants.image_similarity(img1, img2, metric_type=critere, fixed_mask=None, moving_mask=path_mask)  # Calcul similarite selon critere imposé et masque de l'image mouvante si donné


def Recalage_atlas(atlas_fix, img_mouv, type_transfo, interpolator):
    """

    :param atlas_fix: image fixe ou de reference vers laquelle on va recaler 'ici l'atlas)
    :param img_mouv: image à recalé vers l'image referente (ici image sujet)
    :param type_transfo: type de fonction de recalage (rigide, affine, SyN etc)
    :param interpolator: type d'interpolation lors du recalage (proche voisins, lineaires etc)
    :return:
    """
    warp_sub = ants.registration(atlas_fix, img_mouv, type_of_transform=type_transfo)  #
    return ants.apply_transforms(atlas_fix, img_mouv, transformlist=warp_sub['fwdtransforms'], interpolator=interpolator)


def SAVE_Transfo_rec_mat(atlas_fix, img_mouv, type_transfo, file_transfo_direct, file_transfo_inv, name_sujet, name_atlas):
    """

    :param atlas_fix:
    :param img_mouv:
    :param type_transfo:
    :param file_transfo_direct:
    :param file_transfo_inv:
    :param name_sujet:
    :param name_atlas:
    :return:
    """
    path_file_transfo_direct = creation_chemin_fichier_mat(file_transfo_direct, name_sujet, name_atlas)
    path_file_transfo_inv = creation_chemin_fichier_mat(file_transfo_inv, name_sujet, name_atlas)
    ants.registration(atlas_fix, img_mouv, type_of_transform=type_transfo, outprefix=path_file_transfo_direct + '_direct_')
    ants.registration(img_mouv, atlas_fix, type_of_transform=type_transfo, outprefix=path_file_transfo_inv + '_Inverse_')
    return path_file_transfo_direct, path_file_transfo_inv


def path_abs_sujet_to_fichier_repertorie_sujet(tab_path):
    """

    :param tab_path:
    :return:
    """
    repertoire = [os.path.dirname(path) for path in tab_path]
    fichier = [os.path.basename(path) for path in tab_path]
    return repertoire, fichier


def Enregistrer_img_ants_en_nifit(img, path_repertoire, nom_img):
    """

    :param img:
    :param path_repertoire:
    :param nom_img:
    :return:
    """
    ants.image_write(img, os.path.join(path_repertoire, nom_img))


def tab2d_atlas_sim_critere(lignes_atlas,criteres):
    """

    :param lignes_atlas:
    :param criteres:
    :return:
    """
    tab2D = np.zeros((len(lignes_atlas), 3), dtype=object)
    tab2D[:, 0] = lignes_atlas
    tab2D[:, 2] = np.array(criteres[0])
    return tab2D


def recupAtlas_to_tableau_simil(lignes_atlas, criteres, path_atlas, sujet, sujet_repertoire, type_transfo, interpolation, mask = None):
    """

    :param lignes_atlas:
    :param criteres:
    :param path_atlas:
    :param sujet:
    :param sujet_repertoire:
    :param type_transfo:
    :param interpolation:
    :param mask:
    :return:
    """
    tab2D = tab2d_atlas_sim_critere(lignes_atlas, criteres)
    sujet_ants = ants.image_read((os.path.join(sujet_repertoire, sujet)))
    for i in range(len(tab2D[:, 0])):
        Atlas_recherche = ants.image_read((os.path.join(path_atlas, tab2D[i, 0])))
        Sujet_Warped = Recalage_atlas(Atlas_recherche, sujet_ants, type_transfo, interpolation)
        for critere in criteres:
            similarity = calcul_similarity_ants(Atlas_recherche, Sujet_Warped, critere, mask)
            tab2D[i, 1] = similarity
    plot_sujet_by_atlas_simil(tab2D[:, 0], tab2D[:, 1], sujet)
    return tab2D


def atlas_du_bon_age(lignes_atlas, criteres, path_atlas, sujet, sujet_repertoire, type_transfo, interpolation, mask = None):
    """

    :param lignes_atlas:
    :param criteres:
    :param path_atlas:
    :param sujet:
    :param sujet_repertoire:
    :param type_transfo:
    :param interpolation:
    :param mask:
    :return:
    """
    tab_similarity = recupAtlas_to_tableau_simil(lignes_atlas, criteres, path_atlas, sujet, sujet_repertoire, type_transfo, interpolation, mask)
    indice_val_max = np.argmax(np.abs(tab_similarity[:, 1].astype(float)))
    nom_max = tab_similarity[indice_val_max, 0]
    return nom_max


def recup_bon_atlas_avc_transfos(lignes_atlas, criteres, path_atlas, sujet, sujet_repertoire, type_transfo, interpolation, file_transfo_direct, file_transfo_inv, mask = None):
    """

    :param lignes_atlas:
    :param criteres:
    :param path_atlas:
    :param sujet:
    :param sujet_repertoire:
    :param type_transfo:
    :param interpolation:
    :param file_transfo_direct:
    :param file_transfo_inv:
    :param mask:
    :return:
    """
    bon_atlas = atlas_du_bon_age(lignes_atlas, criteres, path_atlas, sujet, sujet_repertoire, type_transfo, interpolation,mask)
    sujet_ants = ants.image_read((os.path.join(sujet_repertoire, sujet)))
    atlas_ants = ants.image_read((os.path.join(path_atlas, bon_atlas)))
    path_trf_direct, path_trf_inv = SAVE_Transfo_rec_mat(atlas_ants, sujet_ants, type_transfo, file_transfo_direct, file_transfo_inv, sujet, bon_atlas)
    return bon_atlas, path_trf_direct, path_trf_inv


def creation_chemin_nom_img(path_repertoire_output, img_name, suffix_nom_image: str):
    """

    :param path_repertoire_output:
    :param img_name:
    :param suffix_nom_image:
    :return:
    """
    nom_initial, fin = (img_name[:-7], ".nii.gz") if img_name.endswith(".nii.gz") else os.path.splitext(img_name)
    return os.path.join(path_repertoire_output, f"{nom_initial}_{suffix_nom_image}")


def creation_chemin_fichier_mat(path_repertoire_output,img_name, atlas_name):
    """

    :param path_repertoire_output:
    :param img_name:
    :param atlas_name:
    :return:
    """
    nom_initial, fin = (img_name[:-7], ".gz") if img_name.endswith(".nii.gz") else os.path.splitext(img_name)
    nom_2, fin2 = (atlas_name[:-7], ".gz") if atlas_name.endswith(".nii.gz") else os.path.splitext(atlas_name)
    return os.path.join(path_repertoire_output, f"{nom_initial}_to_{nom_2}")


def plot_sujet_by_atlas_simil(list1,list2,sujet):
    """

    :param list1:
    :param list2:
    :param sujet:
    :return:
    """
    numero_atlas_x = extraction_numero_atlas(list1)
    similarite_abs = np.abs(list2.astype(float))
    plt.figure(figsize = (10, 6))
    plt.plot(numero_atlas_x, similarite_abs, marker = 'o')
    plt.title(f" Valeurs de similarité pour chaque atlas pour le sujet {sujet[:-7]}", fontsize = 11, pad = 20)
    plt.xlabel('l\'age de l\'atlas en semaine ', fontsize = 12)
    plt.ylabel('valeur d\'information mutuelle', fontsize = 12)
    plt.grid(True)
    plt.tight_layout()
    plt.show()



def extraction_numero_atlas(list_atlas):
    """

    :param list_atlas:
    :return:
    """
    list_num = []
    for atlas in list_atlas:
        nom, fin = (atlas[:-7], ".nii.gz") if atlas.endswith(".nii.gz") else os.path.splitext(atlas)
        numero_atlas = nom.split('STA')[1]
        list_num.append(numero_atlas)
    return list_num



def extraction_numero_sujet(list_sujet):
    """

    :param list_sujet:
    :return:
    """
    list_nums = []
    for sujet in list_sujet :
        num, fin = (sujet[:-49], "_acq-haste_rec-nesvor_desc-aligned_T2w_rot.nii.gz") if sujet.endswith("_acq-haste_rec-nesvor_desc-aligned_T2w_rot.nii.gz") else os.path.splitext(sujet)
        list_nums.append(num)
    return list_nums


def creation_data_frame_sujet_by_best_atlas(list_sujet, list_atlas):
    """

    :param list_sujet:
    :param list_atlas:
    :return:
    """
    age_atlas = extraction_numero_atlas(list_atlas)
    nums_sujet = extraction_numero_sujet(list_sujet)
    reel_age = ["28.4", "20.8","32.3", "29", "24.4"]
    data = {'les numeros du sujet       ': nums_sujet, ' age (semaine) du meilleur atlas    ': age_atlas, 'age estime (pendant acquisition irm)': reel_age}
    df = pd.DataFrame(data)
    pd.set_option('display.colheader_justify', 'center') #On force le centrage des colonnes
    print(df.to_string(index=False))
