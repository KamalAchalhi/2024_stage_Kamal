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

    :param nom_general_sujet:
    :param repertoire_sujet_segm:
    :param pattern_sous_repertoire_by_sujet:
    :return:
    """
    file_paths = []
    pattern = re.compile(nom_general_sujet)

    if repertoire_sujet_segm :
        for file in os.listdir(repertoire_sujet_segm):
             if pattern.match(file):
                    file_paths.append(os.path.join(repertoire_sujet_segm ,file))
    elif pattern_sous_repertoire_by_sujet:
        path_pattern = re.compile(pattern_sous_repertoire_by_sujet)
        for root, _, files in os.walk("/envau/work/meca/users/2024_Kamal/real_data/lastest_nesvor/"):
            if path_pattern.search(root):
                for file in files :
                    if pattern.match(file):
                        file_paths.append(os.path.join(root, file))
    return sorted(file_paths)


def copy_info_geo(path_img_input, path_img_input_copied):
    """

    :param path_img_input:
    :param path_img_input_copied:
    :return:
    """
    img_copied = nib.load(path_img_input_copied)
    print(np.shape(img_copied))
    img_input = nib.load(path_img_input).get_fdata()
    print(np.shape(img_input))
    return nib.Nifti1Image(img_input, img_copied.affine, img_copied.header)


def SWAP_COPY_INFO_SAVE(path_img_input, path_img_rot_nifti):
    """

    :param path_img_input:
    :param path_img_rot_nifti:
    :return:
    """
    img_input = nib.load(path_img_input)
    img_input_array = img_input.get_fdata()
    original_data_type = img_input.get_data_dtype()
    img_input_array_rot = np.transpose(img_input_array, (0, 1, 2))[::1, ::1, ::-1]
    img_input_array_rot = img_input_array_rot.astype(original_data_type)
    img_input_rot_nifti = nib.Nifti1Image(img_input_array_rot, img_input.affine, img_input.header)
    nib.save(img_input_rot_nifti, path_img_rot_nifti)


def creation_PATH_pour_fichier_swaper(path_sujet, repertoire_output):
    """

    :param path_sujet:
    :param repertoire_output:
    :return:
    """
    fichier = os.path.basename(path_sujet)
    print(fichier)
    nom_initial, fin = (fichier[:-7], ".nii.gz") if fichier.endswith(".nii.gz") else os.path.splitext(fichier)
    return os.path.join(repertoire_output, f"{nom_initial}_rot{fin}")


def Parcours_dossier_only_data_match(Path, nom_caracteristic):
    """

    :param Path:
    :param nom_caracteristic:
    :return:
    """
    pattern = re.compile(nom_caracteristic)
    return sorted(f for f in os.listdir(Path) if pattern.match(f))


def calcul_similarity_ants(img1, img2, critere, path_mask = None):
    """

    :param img1:
    :param img2:
    :param critere:
    :param path_mask:
    :return:
    """
    return ants.image_similarity(img1, img2, metric_type=critere, fixed_mask=None, moving_mask=path_mask)


def Recalage_atlas(atlas_fix, img_mouv, type_transfo, interpolator):
    """

    :param atlas_fix:
    :param img_mouv:
    :param type_transfo:
    :param interpolator:
    :return:
    """
    warp_sub = ants.registration(atlas_fix, img_mouv, type_of_transform=type_transfo)
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
