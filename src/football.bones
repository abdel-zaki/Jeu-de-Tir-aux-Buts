# Blender3D Bones File: football.bones
transl ROOT_transl
-1.5105 -2.9005 2.1803
rot_quat Colonne_rot
0.7427 0.6696 -0.0008 0.001
bone Colonne
0.9963
parent_bone NULL

  transl Colonne_transl
  0 0.9963 0
  rot_quat Pectoraux_L_rot
  0.5379 -0.4745 -0.4434 -0.5375
  bone Pectoraux_L
  0.3335
  parent_bone Colonne

    transl Pectoraux_L_transl
    0 0.3335 0
    rot_quat Bras_L_rot
    0.703 -0.4628 0.4573 -0.2874
    bone Bras_L
    0.4527
    parent_bone Pectoraux_L

      transl Bras_L_transl
      0 0.4527 0
      rot_quat Avant_Bras_L_rot
      0.9591 0.1889 0.1132 -0.178
      bone Avant_Bras_L
      0.4692
      parent_bone Bras_L

        transl Avant_Bras_L_transl
        0 0.4692 0
        rot_quat Main_L_rot
        0.9934 0.0157 -0.0856 0.0752
        bone Main_L
        0.1985
        parent_bone Avant_Bras_L

        bone_end Main_L

      bone_end Avant_Bras_L

    bone_end Bras_L

  bone_end Pectoraux_L

  transl Colonne_transl
  0 0.9963 0
  rot_quat Pectoraux_R_rot
  0.539 -0.4732 0.4447 0.5364
  bone Pectoraux_R
  0.3283
  parent_bone Colonne

    transl Pectoraux_R_transl
    0 0.3283 0
    rot_quat Bras_R_rot
    0.6398 -0.435 -0.4365 0.4593
    bone Bras_R
    0.4766
    parent_bone Pectoraux_R

      transl Bras_R_transl
      0 0.4766 0
      rot_quat Avant_Bras_R_rot
      0.9499 0.3096 0.0354 -0.0246
      bone Avant_Bras_R
      0.4814
      parent_bone Bras_R

        transl Avant_Bras_R_transl
        0 0.4814 0
        rot_quat Main_R_rot
        0.9924 -0.0961 -0.0171 -0.0756
        bone Main_R
        0.1628
        parent_bone Avant_Bras_R

        bone_end Main_R

      bone_end Avant_Bras_R

    bone_end Bras_R

  bone_end Pectoraux_R

  transl Colonne_transl
  0 0.9963 0
  rot_quat Tete_rot
  0.0006 0.0112 0.9999 -0.0043
  bone Tete
  0.5009
  parent_bone Colonne

  bone_end Tete

bone_end Colonne

transl ROOT_transl
-1.5105 -2.9005 2.1803
rot_quat Cuisse_L_rot
0.2395 0.0142 0.6664 -0.7059
bone Cuisse_L
0.7077
parent_bone NULL

  transl Cuisse_L_transl
  0 0.7077 0
  rot_quat Tebia_L_rot
  0.9881 -0.0934 0.0818 -0.0904
  bone Tebia_L
  0.7043
  parent_bone Cuisse_L

    transl Tebia_L_transl
    0 0.7043 0
    rot_quat Pied_L_rot
    0.9187 0.3882 0.0717 -0.0106
    bone Pied_L
    0.2918
    parent_bone Tebia_L

    bone_end Pied_L

  bone_end Tebia_L

bone_end Cuisse_L

transl ROOT_transl
-1.5105 -2.9005 2.1803
rot_quat Cuisse_R_rot
0.1295 0.0201 -0.6434 0.7543
bone Cuisse_R
0.6974
parent_bone NULL

  transl Cuisse_R_transl
  0 0.6974 0
  rot_quat Tebia_R_rot
  0.9834 -0.0973 -0.0933 0.1215
  bone Tebia_R
  0.6937
  parent_bone Cuisse_R

    transl Tebia_R_transl
    0 0.6937 0
    rot_quat Pied_R_rot
    0.9218 0.3625 0.1285 -0.0488
    bone Pied_R
    0.2804
    parent_bone Tebia_R

    bone_end Pied_R

  bone_end Tebia_R

bone_end Cuisse_R

