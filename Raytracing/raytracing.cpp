#include "raytracing.h"

Primitives Node::prim;

float clamp01(float v) {
  if (v < 0.0f)
    return 0.0f;
  if (v > 1.0f)
    return 1.0f;
  return v;
}

Vec3 Color(const Vec3 &c) {
  return Vec3(clamp01(c.x), clamp01(c.y), clamp01(c.z));
}

/*
 *  /////////////////////////////////////////
 *  /////////////// PRIMITIVES //////////////
 *  /////////////////////////////////////////
 */
Node::Node(const Mat4 &m, const Vec3 &c, float sp, float tr)
    : transfo(m), inv_transfo(glm::inverse(m)), Col(c), spec(sp), transp(tr) {}

void Cube::draw_gl() const {
  Node::prim.draw_cube(this->transfo, this->Col);
  // DC
}
void Sphere::draw_gl() const {
  Node::prim.draw_sphere(this->transfo, this->Col);
  // DC
}
void Cylinder::draw_gl() const {
  Node::prim.draw_cylinder(this->transfo, this->Col);
  // DC
}

/*
 * Origin : l'origine du rayon
 * Dir : la direction du rayon
 * I : Informations sur le point d'intersection
 * return true en cas d'intersection
 */
bool Cube::intersecte(const Vec3 &Origin, const Vec3 &Dir, Inter *I) {
  // AC

  // Transformer le rayon pour tester l'intersection avec un cube centré en
  // (0,0,0) (utilisation de inv_transfo)

  // Tester l'intersection avec les 6 faces

  // Retransformer le point d'intersection pour le ramener dans la scène
  // (utilisation de transfo)
}

/*
 * P : point sur l'objet auquel on veut calculer la normale
 * retourne la normal normalisé
 */
Vec3 Cube::normal(const Vec3 &P) {
  // AC. Pensez à normaliser la normale une fois le vecteur calculé.

  // Calcul de la normale avec un cube centré en (0,0,0) (utilisation de
  // inv_transfo)

  // Un ramène la normale dans la scène (utilisation de transfo)
}

/*
 * Origin : l'origine du rayon
 * Dir : la direction du rayon
 * I : Informations sur le point d'intersection
 * return true en cas d'intersection
 */
bool Sphere::intersecte(const Vec3 &Origin, const Vec3 &Dir, Inter *I) {
  // AC
  // Meme principe que pour le cube
}

/*
 * P : point sur l'objet auquel on veut calculer la normale
 * retourne la normal normalisé
 */
Vec3 Sphere::normal(const Vec3 &P) {
  // AC
  // Meme principe que pour le cube
}

/*
 * Origin : l'origine du rayon
 * Dir : la direction du rayon
 * I : Informations sur le point d'intersection
 * return true en cas d'intersection
 */
bool Cylinder::intersecte(const Vec3 &Origin, const Vec3 &Dir, Inter *I) {
  // AC
  // Meme principe que pour le cube
}

/*
 * P : point sur l'objet auquel on veut calculer la normale
 * retourne la normal normalisé
 */
Vec3 Cylinder::normal(const Vec3 &P) {
  // AC
  // Meme principe que pour le cube
}

/*
 *  /////////////////////////////////////////
 *  ////////////////// BVH //////////////////
 *  /////////////////////////////////////////
 */

/*
 * Origin : l'origine du rayon
 * Dir : la direction du rayon
 * I : Informations sur le point d'intersection
 */
void BVH::closestIntersection(const Vec3 &Origin, const Vec3 &Dir, Inter *I) {
  // AC
  // On parcours le BVH (structure d'arbre) en testant l'intersection avec le rayon (Origin et Dir en paramètres)
  // avec les boites englobante des noeuds intermédiaires du BVH puis les primitives (Nodes) sur les feuilles
  // On cherche l'intersection LA PLUS PROCHE de l'origine du rayon !
  // On stocke les informations d'intersection dans la struct Inter I
}

/*
 * Origin : l'origine du rayon (de la lumière)
 * Dir : la direction du rayon (de la lumière)
 * sha : coefficient d'ombre (entre 0.0 et 1.0)
 * return true en cas d'intersection
 */
void BVH::intersecteShadow(const Vec3 &Origin, const Vec3 &Dir, float &sha) {
  // AC
  // On parcours le BVH (structure d'arbre) en testant l'intersection avec le rayon (Origin et Dir en paramètres)
  // avec les boites englobante des noeuds intermédiaires du BVH puis les primitives (Nodes) sur les feuilles
  // Ici (contrairement à la fct BVH::intersecte()), on cherche juste si il y a intersection ou pas !

  // Si il y a intersection alors le point d'origine est dans l'ombre car il y a un autre objet entre lui et la source de lumière.

  // ATTENTION tout de même : si l'objet intersecté est transparent, on cumule son coefficient de transparence
  // et on continue le test si celui n'a pas atteint une valeur de 1.0 !

  // On met à jour le coefficient d'ombre sha pris en paramètre si il y a des intersections et selon le coefficient de transparence des objet intersectés.
}

const float Node::Epsilon = 0.0001f;

/*
 *  /////////////////////////////////////////
 *  /////////////// RAY TRACER //////////////
 *  /////////////////////////////////////////
 */

RTracer::RTracer(SimpleViewer *s) : sv(s), depth(0) {}

void RTracer::add_cube_bvh(BVH *b, const Mat4 &m, const Vec3 &color, float spec,
                           float tr) {
  b->nodes.push_back(new Cube(m, color, spec, tr));
  // DC
}
void RTracer::add_sphere_bvh(BVH *b, const Mat4 &m, const Vec3 &color,
                             float spec, float tr) {
  b->nodes.push_back(new Sphere(m, color, spec, tr));
  // DC
}
void RTracer::add_cylinder_bvh(BVH *b, const Mat4 &m, const Vec3 &color,
                               float spec, float tr) {
  b->nodes.push_back(new Cylinder(m, color, spec, tr));
  // DC
}
void RTracer::add_sponge_bvh(BVH *b, const Mat4 &m, const Vec3 &color,
                               float spec, float tr, int r) {
    float multiplier = 2./3.+0.1;
    float taille = 1./3.;
    if (r == 0) {
        //draw cube and sphere
        add_cube_bvh(b, m, color, spec, tr);
        return;
    }
    for (float i=-multiplier; i<=multiplier;i+=multiplier)
        for (float j=-multiplier; j<=multiplier; j+=multiplier*2)
            for (float k=-multiplier; k<=multiplier; k+=multiplier*2)
                add_sponge_bvh(b, m*translate(i, j, k)*scale(taille), color, spec, tr, r-1);

    for (float i=-multiplier; i<=multiplier; i+=multiplier*2) {
        add_sponge_bvh(b, m*translate(i, multiplier, 0)*scale(taille), color, spec, tr, r-1);
        add_sponge_bvh(b, m*translate(i, -multiplier, 0)*scale(taille), color, spec, tr, r-1);
        add_sponge_bvh(b, m*translate(i, 0, multiplier)*scale(taille), color, spec, tr, r-1);
        add_sponge_bvh(b, m*translate(i, 0, -multiplier)*scale(taille), color, spec, tr, r-1);
    }
    multiplier*=3;
    for(float i=-multiplier; i<=multiplier; i+=multiplier*2) {
        add_sphere_bvh(b, m*scale(taille)*translate(i, 0, 0), color, spec, tr);
        add_sphere_bvh(b, m*scale(taille)*translate(0, i, 0), color, spec, tr);
        add_sphere_bvh(b, m*scale(taille)*translate(0, 0, i), color, spec, tr);
    }
}

Vec3 RTracer::ColorRayBVH(const Vec3 &Origin, const Vec3 &Dir, int rec) {
  Inter I;
  I.a_min = std::numeric_limits<float>::max();
  I.node = nullptr;

  bvh->closestIntersection(Origin, Dir, &I);

  if (I.node != nullptr) {
    const Vec3 &D = I.Dir;

    // Position
    Vec3 Po = I.Pos;

    // Light direction
    Vec3 Lu = Vec3(0, 0, 0);
    // AC : calculer Lu, la direction de la lumière

    // Normal
    Vec3 No = I.node->normal(Po);

    // Shadow
    float sha = 0.0f;
    // AC : calculer le coefficient d'ombre

    // Lambert (diffuse) BRDF
    float lambert = 0.0f;
    // AC : calculer la diffusion de l'objet avec -> (1 - sha) * col ∗ coeff ∗ max(0, No · Lu)

    // ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // La suite est déjà complétée, il sagit d'autres comportements de la lumière et des objets : spécularité, transparence et réflexions
    // ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // View direction for specular BRDF
    Vec3 inci = glm::normalize(D);
    Vec3 R = inci - 2.0f * glm::dot(No, inci) * No;
    Vec3 RR = glm::normalize(R);

    if (rec > 0) {
      Vec3 col2(0);
      if (I.node->spec > 0)
        col2 = ColorRayBVH(Po, RR, rec - 1);

      Vec3 col3(0);
      if (I.node->transp > 0)
        col3 = ColorRayBVH(Po, D, rec - 1);

      Vec3 final = (1.0f - I.node->transp) *
                       ((1.0f - I.node->spec) * I.node->Col * lambert +
                        I.node->spec * col2) +
                   I.node->transp * col3;

      // Specular BRDF
      final += (1.0f - sha) * I.node->spec *
               Vec3(std::pow(std::max(0.0f, glm::dot(RR, Lu)), 40.0f));

      return Color(final);
    }

    Vec3 final = I.node->Col * lambert;

    // Specular BRDF
    final += (1.0f - sha) * I.node->spec *
             Vec3(std::pow(std::max(0.0f, glm::dot(RR, Lu)), 40.0f));

    return Color(final);
  }
  return Vec3(0, 0, 0);
}

QLabel *RTracer::CalcImage(int depth) {
  std::cout << "Calcul lancer de rayon avec " << depth << " rebonds -- sur "
            << this->nbthr << " threads ..." << std::endl;

  QImage im(this->sv->width(), this->sv->height(), QImage::Format_RGB32);
  Vec4 SHADER_POSLUM(300, 1000, 1000, 1);
  this->posLum =
      Vec3(glm::inverse(this->sv->getCurrentModelViewMatrix()) * SHADER_POSLUM);
  this->depth = depth;

  int x, y;
  auto start = std::chrono::system_clock::now();
#pragma omp parallel for private(x) schedule(dynamic) num_threads(this->nbthr)
  for (y = 0; y < this->sv->height(); ++y) {
    for (x = 0; x < this->sv->width(); ++x) {
      qglviewer::Vec Pq = this->sv->camera()->unprojectedCoordinatesOf(
          qglviewer::Vec(x, y, 0.0));
      qglviewer::Vec Qq = this->sv->camera()->unprojectedCoordinatesOf(
          qglviewer::Vec(x, y, 1.0));
      Vec3 P(Pq[0], Pq[1], Pq[2]);
      Vec3 D(Qq[0] - Pq[0], Qq[1] - Pq[1], Qq[2] - Pq[2]);
      Vec3 C = this->ColorRayBVH(P, D, this->depth);
      // ^ ou ColorRay(...) : lance un rayon d'origine P et de direction D,
      // rebondissant this->depth fois.
      im.setPixel(
          x, y,
          QColor(int(C.r * 255.0f), int(C.g * 255.0f), int(C.b * 255.0f))
              .rgb());
    }
  }
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed = end - start;
  std::cout << "Elapsed time: " << elapsed.count() << "s" << std::endl;

  QLabel *label_img = new QLabel(this->sv);
  label_img->setWindowFlags(Qt::Window);
  label_img->setPixmap(QPixmap::fromImage(im));
  return label_img;
  // DC
}

void draw_prim_bvh(BVH *b) {
  for (const auto *n : b->nodes)
    n->draw_gl();
  for (auto *c : b->children)
    draw_prim_bvh(c);
  // DC
}
