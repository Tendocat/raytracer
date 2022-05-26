#ifndef __RAYTRACING_H__
#define __RAYTRACING_H__

#include "QGLViewer/simple_viewer.h"
#include "matrices.h"
#include "primitives.h"
#include <thread>

#include <QApplication>
#include <QLabel>

const Vec3 ROUGE = {1, 0, 0};      // 0
const Vec3 VERT = {0, 1, 0};       // 1
const Vec3 BLEU = {0, 0, 1};       // 2
const Vec3 JAUNE = {1, 1, 0};      // 3
const Vec3 CYAN = {0, 1, 1};       // 4
const Vec3 MAGENTA = {1, 0, 1};    // 5
const Vec3 BLANC = {1, 1, 1};      // 6
const Vec3 GRIS = {0.5, 0.5, 0.5}; // 7
const Vec3 NOIR = {0, 0, 0};       // 8
const Vec3 ROUGE2 = {0.5, 0, 0};   // 9
const Vec3 VERT2 = {0, 0.5, 0};    // 10
const Vec3 ORANGE = {1, 0.5, 0};   // 11

class Node;
struct Inter {
  float a_min; // distance de l'intersection par rapport à l'origine du rayon.
  Node *node;  // noeud intercepté, if any.
  Vec3 Pos; // point d'impact du rayon = o + dt (avec o : l'origine du rayon; d : sa direction; et t : la distance de l'intersection par rapport à l'origine du rayon)
  Vec3 Dir; // direction du rayon. (ps : ctrl+shift+R sur qtcreator pour changer
            // le nom des variables).
};

class Node {
public:
  static const float Epsilon;
  Mat4 transfo; // matrice de transformation.
  Mat4 inv_transfo;
  Vec3 Col;     // couleur.
  float spec;   // indice de spécularité.
  float transp; // indice de transparence.
  Node(const Mat4 &m, const Vec3 &c, float sp, float tr);
  inline virtual ~Node() {}
  virtual void draw_gl() const = 0; // la fonction permettant de dessiner la
                                    // primitive dans la vue temps réel.
  virtual bool intersecte(const Vec3 &P, const Vec3 &V,
                          Inter *I) = 0;  // test d'intersection de la primitive
                                          // (modifier contenu de I).
  virtual Vec3 normal(const Vec3 &P) = 0; // la normale à la primitive

  inline bool isEqual(float a, float b) {
    return fabs(a-b) < Epsilon;
  }

  static Primitives prim;
};

class Cube : public Node {
public:
  inline Cube(const Mat4 &_transfo, const Vec3 &_color, float _spec,
              float _transp)
      : Node(_transfo, _color, _spec, _transp) {}
  void draw_gl() const;
  bool intersecte(const Vec3 &Origin, const Vec3 &Dir, Inter *I);
  Vec3 normal(const Vec3 &P);
};

class Sphere : public Node {
public:
  inline Sphere(const Mat4 &_transfo, const Vec3 &_color, float _spec,
                float _transp)
      : Node(_transfo, _color, _spec, _transp) {}
  void draw_gl() const;
  bool intersecte(const Vec3 &Origin, const Vec3 &Dir, Inter *I);
  Vec3 normal(const Vec3 &P);
};

class Cylinder : public Node {
public:
  inline Cylinder(const Mat4 &_transfo, const Vec3 &_color, float _spec,
                  float _transp)
      : Node(_transfo, _color, _spec, _transp) {}
  void draw_gl() const;
  bool intersecte(const Vec3 &Origin, const Vec3 &Dir, Inter *I);
  Vec3 normal(const Vec3 &P);
};

class BVH {
public:
  Mat4 transfo;
  Mat4 inv_transfo;
  std::vector<BVH *> children;
  std::vector<Node *> nodes;

  inline BVH(const Mat4 &m) : transfo(m), inv_transfo(glm::inverse(m)) {}

  inline BVH *add_child(const Mat4 &m) {
    BVH *b = new BVH(m);
    children.push_back(b);
    return b;
  }

  inline void add_node(Node *n) { nodes.push_back(n); }

  // Parcours du BVH pour calculer l'intersection entre les rayons primaires (rayon qui passe par un pixel et se dirige vers la scène)
  // et les primitives de la scène.
  void closestIntersection(const Vec3 &P, const Vec3 &V, Inter *I);

  // Parcours du BVH pour calculer les ombres. Ici le rayon attendu à pour origine un point d'intersection
  // d'un rayon primaire sur une primitive, et une direction vers la source lumineuse
  void intersecteShadow(const Vec3 &Origin, const Vec3 &Dir, float &sha);
  static float norme(Vec3 a, Vec3 b);
};

class RTracer {
public:
  BVH *bvh;
  SimpleViewer *sv;
  std::vector<Node *> primitives;
  int depth;
  Vec3 posLum;  // Position de la source lumineuse
  static const int nbthr = 8;

  RTracer(SimpleViewer *s);

  void add_cube_bvh(BVH *b, const Mat4 &m, const Vec3 &color, float spec,
                    float tr);
  void add_sphere_bvh(BVH *b, const Mat4 &m, const Vec3 &color, float spec,
                      float tr);
  void add_cylinder_bvh(BVH *b, const Mat4 &m, const Vec3 &color, float spec,
                        float tr);
  void add_sponge_bvh(BVH *b, const Mat4 &m, const Vec3 &color, float spec,
                      float tr, int r);
  void add_apoll_bvh(BVH *b, const Mat4 &m, const Vec3 &color, int r);

  QLabel *CalcImage(int depth); // calcule l'image par lancer de rayons. Renvoie
                                // de quoi l'afficher.

  // v   à utiliser si vous voulez tester un rendu de la scène avant
  // d'implémenter un bvh.
  Vec3 ColorRay(const Vec3 &P, const Vec3 &V, int depth);
  Vec3 ColorRayBVH(const Vec3 &Origin, const Vec3 &Dir, int rec);
};

void draw_prim_bvh(BVH *b);

#endif
