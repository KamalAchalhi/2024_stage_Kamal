# code développé par Kamal durant son stage en 2024

## installation
### cloner le repos

### installer un environnement virtuel
module load all anaconda

conda create --name stage_Kamal python=3.8

conda activate stage_Kamal

pip install antspyx

### organisation du repos
#### Dossier CODE
#####tools : Ce fichier python contient toutes les fonctions que j'utilise lors des diffèrentes étapes de la pipeline :
Fonction de récuperéation de fichier selon pattern.
fonction de création de chemin (récupère le path d'entrée et  modifie le nom du fichier pour qu'il puisse être adapté au processus qui a été utilisé (ajout de ROT pour fichier après 'swap', 
fusionnée des noms apres recalage et indiquer direction de recalage ex:  "sub-0002_ses-0002_acq-haste_rec-nesvor_desc-aligned_T2w_rot_STA21_all_reg_LR_dilM.nii.gz")).
Fonction qui recale image sujet à une liste donnée d'atlas et qui estime la similarité.
Fonction qui cherche et retourne l'atlas qui maximise la similarité.
Fonction qui enregistre la transformé inverse de recalage dans un output choisis, retourne le path.
Fonction qui 'swap' et copie et enregistre les informations géometriques.


#####script1 :
ce script permet de recuperer les images de nos sujets à partir d'un nom pattern (dont les numeros varient) qu'on lui donne en entrée,
d'un chemin pattern (dont le numéro varie selon le sujet). 
Je crée des chemins pour les images des sujets après transposition/basculement ("swap").
Je "swap" et copies les informations géométriques de chacune de mes images sujets recupérés.
Pour chacun de mes sujets, j'utilise la "fonction recalage par les atlas" de mon tools, et je recupère en sortir dans des listes, l'atlas du bon âge, le path de chaque transformation inverse de recalage.

#####script2 :
Ce script charge les sorties du 1er script, il récupère les atlas binarisés correspondant aux atlas déterminés prcédemment. 
Sur une boucle, il applique les transformations inverses obtenues dans le script1 du sujet vers l'atlas, obtenant donc le sens inverse, 
et l'applique sur les atlas binaires vers le sujet correspondant. Il crée ensuite un chemin adaptés vers une sortie et un nom d'image adaptés

#####scriptFIN : 
Ce script récupère les chemins vers les atlas binaire recalés dans l'espace sujet ainsi que les chemins vers les sujets segmentés.
un swaping est réalisés sur les sujets segmentés avant d'aditionner chaques images chargés en tableau numpy aux atlas d'hémisphères gauche droite recalés dans script2. 
Les parties qui dépassent sont effacés, le format des niveau de gris sont imposés de manière à être entière et les tableaus numpy qui contiennent la combinaison sont sauvé en format nifiti.


#### Le dossier variables (contient des listes de chemins de fichiers sauvés) : 
Ceux des sujets anatomiques que l'on a transposés et inversés selon certains axes (swaping).

Ceux des atlas d'âge adaptés pour chaque sujets mis en entrée dans le script1.

Ceux des transformations direct et inverse obtenue dans le scripts1 entre chaque sujet vers son atlas d'âge adaptés.

Ceux des atlas binarisés projetés vers l'espace sujet par le script2 en utilisant les transformations inverse récupérées.



### Tutoriel
Pour appliquer le pipeline à une image, suivre la procédure suivante:
Il faut preparer pour le script1 les diffèrents chemins nécessaires pour la récupérations des données. 
On lui donne donc :

Path du repertoire contenant les atlas.

Path du repertoire contenant les images sujets.

Pattern du nom des images sujets.

Pattern du nom des Atlas.

Path de sortie pour les données retournées.

Path de sortie pour stockages des variables réutilisées dans les prochains scripts.

Une fois le script1 excecuté, pour initialiser le script2, il faut donner:

Path des variables (pour récupérer les données retournées par le script1).

Path output du script1 qui contient les images de sujet "swapé".

Path du repertoire contenant les Atlas binarisés.

Path de sortie pour stocker les images après le traitement par recalage inverse.

nom des fichiers contenant les diffèrentes variables.

Le scriptFIN nécessite :

Path des variables et nom du fichiers des atlas recalés enregistrés dans le script2.

Path des variables et nom des fichiers contenants la listes des paths des masques des hemisphères gauches, et de celui contenant la listes des sujets.

Path du répertoire des images segmentés des sujets.

Pattern du nom d'une image segmentées.

Path de sortie des images après ce dernier traitement.

![recapitulatif] (/envau/work/meca/users/2024_Kamal/2024_stage_Kamal/thumbnail_SHEMAS_REALISATION.png) 

                        