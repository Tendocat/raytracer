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

    Vec3 o = Vec3(this->inv_transfo*Vec4(Origin, 1));
    Vec3 d = Vec3(this->inv_transfo*Vec4(Dir, 0));

    float ax = (-1-o.x)/d.x;
    float ay = (-1-o.y)/d.y;
    float az = (-1-o.z)/d.z;
    float bx = ( 1-o.x)/d.x;
    float by = ( 1-o.y)/d.y;
    float bz = ( 1-o.z)/d.z;
    float tmp;

    if (ax > bx) {
        tmp = bx;
        bx = ax;
        ax = tmp;
    }
    if (ay > by) {
        tmp = ay;
        ay = by;
        by = tmp;
    }
    if (ax>by || ay>bx)
        return false;
    if (ay > ax)
        ax = ay;
    if (by < bx)
        bx = by;
    if (az > bz) {
        tmp = az;
        az = bz;
        bz = tmp;
    }
    if (ax>bz || az>bx)
        return false;
    if (az > ax)
        ax = az;
    if (bz < bx)
        bx = bz;

    I->Pos = Vec3(this->transfo*Vec4(ax, ay, az, 1));
    I->Dir = Dir;
    I->node = this;
    I->a_min = BVH::norme(I->Pos, Origin);
    return true;
}

/*
 * P : point sur l'objet auquel on veut calculer la normale
 * retourne la normal normalisé
 */
Vec3 Cube::normal(const Vec3 &p) {
  // AC. Pensez à normaliser la normale une fois le vecteur calculé.
    Vec3 r(0,0,0);
    Vec3 P = normalize(Vec3(this->inv_transfo*Vec4(p, 0)));
    if (isEqual(P.x, 1))
        r.x = 1;
    else if (isEqual(P.x, -1))
        r.x = -1;
    else if (isEqual(P.y, 1))
        r.y = 1;
    else if (isEqual(P.y, -1))
        r.y = -1;
    else if (isEqual(P.z, 1))
        r.z = 1;
    else if (isEqual(P.z, -1))
        r.z = -1;
    return normalize(Vec3(this->transfo*Vec4(r,1)));
}

/*
 * Origin : l'origine du rayon
 * Dir : la direction du rayon
 * I : Informations sur le point d'intersection
 * return true en cas d'intersection
 */
bool Sphere::intersecte(const Vec3 &O, const Vec3 &D, Inter *I) {
  // AC
    Vec3 o = Vec3(this->inv_transfo*Vec4(O, 1));
    Vec3 d = Vec3(this->inv_transfo*Vec4(D, 0));
    float a = glm::dot(d, d);
    float b = 2*(glm::dot(o, d));
    float c = glm::dot(o, o) - 1.f;
    float delta = b*b - 4.f*(a*c);
    float r, r1;
    if (delta < 0)
        return false;
    else if (delta < Epsilon) {
        r = -(b)/(2*a);
    } else {
        r1 = ((-b)+sqrt(delta))/(2*a);
        r = ((-b)-sqrt(delta))/(2*a);
        if (r < -Epsilon || r1 < r) {
            r = r1;
        }
        if (r < Epsilon) return false;
    }
    I->Pos = Vec3(this->transfo*Vec4(d.x*r+o.x, d.y*r+o.y, d.z*r+o.z, 1));
    I->Dir = D;
    I->node = this;
    I->a_min = BVH::norme(I->Pos, O);
    return true;
}

/*
 * P : point sur l'objet auquel on veut calculer la normale
 * retourne la normal normalisé
 */
Vec3 Sphere::normal(const Vec3 &P) {
  // AC
  // Meme principe que pour le cube
    return normalize(Vec3(inv_transfo * Vec4(P,1)));
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
    return false;
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
float BVH::norme(Vec3 a, Vec3 b) {
    return glm::length(a - b);
}
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
    Inter test, tmp;

    tmp.a_min = std::numeric_limits<float>::max();
    tmp.node = nullptr;

    if (!Cube(this->transfo, BLANC, 0, 0).intersecte(Origin, Dir, &test))
    return;

    for (auto b : this->children)
        b->closestIntersection(Origin, Dir, I);

    for (auto node : this->nodes) {
        if (node->intersecte(Origin, Dir, &tmp)
        && tmp.a_min < I->a_min) {
            I->Dir = tmp.Dir;
            I->Pos = tmp.Pos;
            I->node = tmp.node;
            I->a_min = tmp.a_min;
        }
    }
}

/*
 * Origin : l'origine du rayon (de la lumière)
 * Dir : la direction du rayon (de la lumière)
 * sha : coefficient d'ombre (entre 0.0 et 1.0)
 */
void BVH::intersecteShadow(const Vec3 &Origin, const Vec3 &Dir, float &sha) {
    Inter test, I;
    Vec3 O = Origin;

    do {
        I.node = nullptr;
        closestIntersection(O, Dir, &I);
        if (I.node != nullptr) {
            sha += 1.0 - I.node->transp;
            O = Vec3(Vec4(I.Pos, 1)*translate(Dir*0.001f));
        }
    } while(I.node != nullptr && sha < 1.0f);

    sha = clamp01(sha);
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

void RTracer::add_sponge_bvh(BVH *B, const Mat4 &m, const Vec3 &color,
                               float spec, float tr, int r) {
    float multiplier = 2./3.+0.1;
    float taille = 1./3.;
    B->add_child(m*scale(2));
    BVH* b = B->children.back();

    if (r == 0) {
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
        add_sphere_bvh(b, m*scale(taille)*translate(i, 0, 0), color, 1, 0);
        add_sphere_bvh(b, m*scale(taille)*translate(0, i, 0), color, 1, 0);
        add_sphere_bvh(b, m*scale(taille)*translate(0, 0, i), color, 1, 0);
    }
}

void RTracer::add_apoll_bvh(BVH *B, const Mat4 &m, const Vec3 &color, int r) {

    float multiplier = 2;
    float taille = 1./2.;
    B->add_child(m*scale(2));
    BVH* b = B->children.back();

    if (r == 0) {
        add_sphere_bvh(b, m, color, 0, 0.3);
        return;
    }
    add_apoll_bvh(b, m*scale(taille)*translate(multiplier, 0, 0), BLEU, r-1);
    add_apoll_bvh(b, m*scale(taille)*translate(0, multiplier, 0), ROUGE, r-1);
    add_apoll_bvh(b, m*scale(taille)*translate(0, 0, multiplier), ROUGE2, r-1);
    add_apoll_bvh(b, m*scale(taille)*translate(-multiplier, 0, 0), VERT, r-1);
    add_apoll_bvh(b, m*scale(taille)*translate(0, -multiplier, 0), VERT2, r-1);
    add_apoll_bvh(b, m*scale(taille)*translate(0, 0, -multiplier), JAUNE, r-1);
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
    Vec3 Lu = normalize(this->posLum - Po);
    // AC : calculer Lu, la direction de la lumière

    // Normal
    Vec3 No = I.node->normal(Po);

    // Shadow
    float sha = 0.0f;
    // AC : calculer le coefficient d'ombre
    bvh->intersecteShadow(Po, Lu, sha);

    float coeff = 1.f;
    // Lambert (diffuse) BRDF
    float lambert = (1-sha) * coeff * fmax(0, glm::dot(Lu, No));

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
  return GRIS;
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
