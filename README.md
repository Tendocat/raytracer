# Geo3D2021

**Antoine Dumoulin**


## travail demandé

### Eponge de menger

La fonction `add_sponge_bvh` permet d'ajouter une éponge à la scène. Elle se comporte
comme pour les primitive, mais prend en plus un paramètre pour le nombre de récursion.
A 0, un seul cube sera créé.

### Ray tracing

Le BVH est fonctionnel, il y a un gain de performance.
J'ai implémenté l'intersection et la normale d'une sphère et d'un cube.

Une matrice de sphère présente l'avancé de mon projet pour le rendu qui se limite au sphères.
Avec 0 rebond, la transparence ne sera pas affecté alors qu'avec 1 ou plus, elle sera apparente.


![sphères](./sphères.png)

## Nouvelle fractale

Je me suis essentiellement concentré sur les sphères. Et c'est pourquoi j'ai choisi une
fractale constitué de sphères. Je me suis inspiré de la baderne d'Apollonius et du
triangle de Sierpiński.

## Scène complète

Les deux grilles du dessus permettent de visualiser les sphères et cubes séparémment.

![scène](./scène.png)
