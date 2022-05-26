
#include <matrices.h>
#include <primitives.h>

#include <QGLViewer/simple_viewer.h>
#include "raytracing.h"

//Légende :
//AC : à compléter
//DC : déjà complet (ne pas toucher sauf en cas d'extrême urgence)
//TC : théoriquement complet (mais modifications possibles en fonction de votre projet)

int main(int argc, char *argv[])
{
	QApplication a(argc, argv);
	SimpleViewer::init_gl();
	SimpleViewer viewer(NOIR,10);
    viewer.setGeometry(10, 10, 800, 800);

	RTracer rt(&viewer);

	// GL init
	viewer.f_init = [&] ()
	{
        Node::prim.gl_init();
        rt.bvh = new BVH(scale(1000));
        float MAX = 5;
        Mat4 grille = translate(0, 5, 0)*rotateX(90)* scale(0.5);
        for (float tr=MAX; tr>0; tr-=1) {
            for (float spec=0; spec<MAX+0.01f; spec+=1) {
                rt.add_sphere_bvh(rt.bvh, grille*translate(tr*2-MAX,spec*2-MAX,0), ROUGE, spec/MAX, tr/MAX);
            }
        }
        rt.add_apoll_bvh(rt.bvh, translate(-4,0,2)*scale(2.5), ROUGE, 2);
        rt.add_sphere_bvh(rt.bvh, translate(0, -6,0)*scale(5, 1, 5), BLANC, 1, 0.4);

        // you can increase the recursion which is the last arg but ray cast will be very long
        rt.add_sponge_bvh(rt.bvh, translate(4,0,0)*scale(3), BLEU, 0, 0, 2);
    };



	viewer.f_draw = [&] ()
	{
		Node::prim.set_matrices(viewer.getCurrentModelViewMatrix(), viewer.getCurrentProjectionMatrix());
        draw_prim_bvh(rt.bvh); //segfault à moins de rajouter des items au bvh
		//TC : le bvh est censé englober toute la scène.
	};


	viewer.f_keyPress = [&] (int key, Qt::KeyboardModifiers /*mod*/)
	{
		//TC : entre 0 et 4, le nombre de rebonds pour votre rendu
		switch(key)
		{
			case Qt::Key_0:
				rt.CalcImage(0)->show();
				break;
			case Qt::Key_1:
				rt.CalcImage(1)->show();
				break;
			case Qt::Key_2:
				rt.CalcImage(2)->show();
				break;
			case Qt::Key_3:
				rt.CalcImage(3)->show();
				break;
			case Qt::Key_4:
				rt.CalcImage(4)->show();
				break;
			default:
				break;
		}
	};


	viewer.show();
	return a.exec();
	//TC
}
