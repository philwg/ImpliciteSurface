#include <iostream>
#include <cmath>
#include <vector>
#include <GL/glut.h>

using namespace std;


/** Declaration des variables
/*
*/
 
const int NMAX  {100};		// Nombre maximum de points saisis par l'utilisateur
const int maxPt {200};		// Nombre de divisions de la zone d'affichage (verticalement et horizontalement)

int width 	 {1200};		// Largeur de la fenêtre OpenGL
int height	 {1200};		// Hauteur de la fenêtre OpenGL
int xStep    {10};			// Le pas courant de discrétisation horizontal 
int yStep    {10};			// Le pas courant de discrétisation vertical
int N 		 {0};			// Le nombre courant de points saisis par l'utilisateur
int selected {-1};			// L'indice du point selectionné courant (-1 si aucun)

int mix_mode {0};			// Le selecteur de mode de mélange des fonctions de potentiel
int pot_mode {1};			// Le selecteur de fonction de potentiel 
int dis_mode {1};			// Le selecteur de mode d'affichage

float seuil  {0.5f};		// Le seuil de potentiel (iso-valeur)
float rayon  {144.0f};		// Le rayon d'influence des fonctions de potentiel

bool droite  {false};		// Le suivi du bouton souris de droite
bool gauche  {false};		// Le suivi du bouton souris de gauche


/** Déclaration d'une structure de point
/*
*/

struct Point {
	
	float x;		// Abscisse du Point
	float y;		// Ordonnée du Point
	int pos;		// Attribut de position pour les sommets de carrés

	Point(float a=0, float b=0)		// Constructeur
	{
		set(a, b);
		this->pos = -1;
	}

	~Point()	// Destructeur
	{}

	/* Méthode d'initialisation / mise à jour des coordonnées */
	void set(float a, float b)
	{
		this->x = a;
		this->y = b;
	}

	/* Méthode de mise à jour de l'attribut de position */
	void setPos(int new_pos)
	{
		this->pos = new_pos;
	}

	/* Méthode de calcul de la distance euclidienne à un autre point */
	float distTo(Point P)
	{
		return sqrt(pow((this->x)-P.x, 2) + pow((this->y)-P.y, 2));
	}

	/* Méthode renvoyant un Point avec des coordonnées modifiées incrémentalement */
	Point add(float delta_x, float delta_y)
	{
		Point P;
		P.set(this->x+delta_x, this->y+delta_y);
		P.setPos(-1);
		return P;
	}

	/* Méthode de renvoi du vertex OpenGL correspondant au Point */
	float* getVertex()
	{
		return new float[2] {this->x, this->y};
	}
};


/** Déclarations des variables à base de Points
/*
*/

Point Pts[NMAX];					// Le tableau des Points saisis par l'utilisateur
vector<vector<Point>> VtoShow;		// Un vecteur de vecteur de Points pour l'affichage


/* Méthodes de reshape du GLUT */
void main_reshape(int w,  int h)
{
	glViewport(0, 0, w, h); 			
	width = w;
	height = h;
	xStep = width/maxPt;		// recalcul du pas de discrétisation horizontal
	yStep = height/maxPt;		// recalcul du pas de discrétisation vertical
	glMatrixMode(GL_PROJECTION); 
	glLoadIdentity(); 
	glOrtho(0.0f, width, 0.0f, height, -1.0f, 1.0f); 
	glMatrixMode(GL_MODELVIEW); 
	glLoadIdentity();
}

/* La méthode de calcul de la fonction de potentiel */
float Threshold(float r) {

	// Murakami ...
	if (pot_mode==1) {
		if (r<rayon) return pow((1-pow(r/rayon, 2)), 2);
		else return 0.0f;
	}

	// Nishimura ...	
	else if (pot_mode==2) {
		if (r<=rayon/3) return 1.0f-3.0f*pow(r/rayon, 2);
		else if (r<=rayon) return (3.0f/2.0f)*pow((1.0f-(r/rayon)), 2);	
		else return 0.0f;					
	}

	// Wyvill ...
	else if (pot_mode==3) {
		if (r<=rayon) return (-4.0f*pow(r/rayon, 6)+17.0f*pow(r/rayon, 4)-22.0f*pow(r/rayon, 2))/9.0f+1.0f;
		else return 0.0f;
	}

	// ou Blinn (a=1 & b=1)			
	else return exp(-pow(r/rayon, 2));									
}

/* La méthode de mixage des influences des différents points saisis */
float Potential_at(Point P)
{
	float level = (mix_mode==-1) ? 12.0f : 0.0f;
	for (int i=0; i<N; i++) {
		if (mix_mode==-1) level = min(level, Threshold(Pts[i].distTo(P)));		// Intersection (min)
		else if (mix_mode==1) level = max(level, Threshold(Pts[i].distTo(P)));	// Union (max)
		else level += Threshold(Pts[i].distTo(P));								// Mélange (Somme)
	}
	return level;	
}

/* Méthode qui teste l'appartenance d'un point à l'ensemble dépassant le seuil */
bool is_Inside(Point P)
{
	return (Potential_at(P) > seuil);
}

/* Méthode qui renvoie le coefficient à utiliser sur le bord d'un carré */
float get_Coef(Point P1, Point P2)
{
	if (dis_mode > 1) {					//------ Affichage proportionnel à la valeur de potentiel sur l'arête
		float v1 = Potential_at(P1);
		float v2 = Potential_at(P2);
		float vm = min(v1, v2);
		float VM = max(v1, v2);
		if ((VM<seuil)||(vm>seuil)) return 0.5f;
		else return (VM-seuil)/(VM-vm);
	}
	else return 0.5f;					//------ Affichage au milieu des arêtes
}

/* Méthode d'affichage des carrés à afficher */
void displaySquares() 
{
	// Pour chaque carré dont au moins un sommet est dans la zone > seuil
	for (vector<vector<Point>>::iterator Sq=VtoShow.begin(); Sq!=VtoShow.end(); Sq++) {
		
		// Fixation des paramètres d'affichage
		glPointSize(1.0f);
		glColor3f(0.96f, 0.48f, 0.24f);
		glPolygonMode(GL_FRONT, GL_FILL);

		// Si le carré n'a qu'un sommet actif ****************************************************
		if (Sq->size()==1) {
			Point P = Sq->at(0);

			// et si l'affichage n'est pas brut
			if (dis_mode > 0) {
				vector<Point> Tri;
				Tri.clear();

				// On calcule les sommets du triangle à afficher suivant sa position
				if (P.pos==0) {
					Tri.push_back(P);
					Tri.push_back(P.add(0.0f, yStep*get_Coef(P, P.add(0.0f, yStep))));
					Tri.push_back(P.add(xStep*get_Coef(P, P.add(xStep, 0.0f)), 0.0f));
				}
				if (P.pos==1) {
					Tri.push_back(P);
					Tri.push_back(P.add(xStep*get_Coef(P, P.add(xStep, 0.0f)), 0.0f));
					Tri.push_back(P.add(0.0f, -yStep*get_Coef(P, P.add(0.0f, -yStep))));
				}
				if (P.pos==2) {
					Tri.push_back(P);
					Tri.push_back(P.add(0.0f, -yStep*get_Coef(P, P.add(0.0f, -yStep))));
					Tri.push_back(P.add(-xStep*get_Coef(P, P.add(-xStep, 0.0f)), 0.0f));
				}
				if (P.pos==3) {
					Tri.push_back(P);
					Tri.push_back(P.add(-xStep*get_Coef(P, P.add(-xStep, 0.0f)), 0.0f));
					Tri.push_back(P.add(0.0f, yStep*get_Coef(P, P.add(0.0f, yStep))));
				}

				// On affiche le triangle
				glBegin(GL_TRIANGLES);
					glVertex2fv(Tri.at(0).getVertex());
					glVertex2fv(Tri.at(1).getVertex());
					glVertex2fv(Tri.at(2).getVertex());
				glEnd();
			}
			else {

				// sinon -> Affichage brut : seulement le point actif
				glBegin(GL_POINTS);
					glVertex2fv(P.getVertex());
				glEnd();
			}
		}

		// Si le carré a deux sommets actifs **************************************************** 
		else if (Sq->size()==2) {

			// et si l'affichage n'est pas brut
			if (dis_mode > 0) {
				Point P1 = Sq->at(0);
				Point P2 = Sq->at(1);

				// Si les deux points sont sur un bord
				if (((P1.pos+P2.pos)%2)==1) {
					vector<Point> Rect;
					Rect.clear();

					// On calcule les sommets du rectangle (ou trapèze) à afficher suivant leurs positions
					if ((P1.pos==0)&&(P2.pos==1)) {
						Rect.push_back(P1);
						Rect.push_back(P2);
						Rect.push_back(P2.add(xStep*get_Coef(P2, P2.add(xStep, 0.0f)), 0.0f));
						Rect.push_back(P1.add(xStep*get_Coef(P1, P1.add(xStep, 0.0f)), 0.0f));
					}
					if ((P1.pos==0)&&(P2.pos==3)) {
						Rect.push_back(P1);
						Rect.push_back(P1.add(0.0f, yStep*get_Coef(P1, P1.add(0.0f, yStep))));
						Rect.push_back(P2.add(0.0f, yStep*get_Coef(P2, P2.add(0.0f, yStep))));
						Rect.push_back(P2);
					}
					if (P1.pos==1) {
						Rect.push_back(P1);
						Rect.push_back(P2);
						Rect.push_back(P2.add(0.0f, -yStep*get_Coef(P2, P2.add(0.0f, -yStep))));
						Rect.push_back(P1.add(0.0f, -yStep*get_Coef(P1, P1.add(0.0f, -yStep))));
					}
					if (P1.pos==2) {
						Rect.push_back(P1);
						Rect.push_back(P2);
						Rect.push_back(P2.add(-xStep*get_Coef(P2, P2.add(-xStep, 0.0f)), 0.0f));
						Rect.push_back(P1.add(-xStep*get_Coef(P1, P1.add(-xStep, 0.0f)), 0.0f));
					}

					// On affiche le rectangle (trapèze en affichage proportionnel)
					glBegin(GL_QUADS);
						glVertex2fv(Rect.at(0).getVertex());
						glVertex2fv(Rect.at(1).getVertex());
						glVertex2fv(Rect.at(2).getVertex());
						glVertex2fv(Rect.at(3).getVertex());
					glEnd();
				}

				// sinon ils sont sur un diagonale
				else {
					vector<Point> Hexa;
					Hexa.clear();

					// On calcule les sommets de l'hexagone à afficher suivant le cas
					if (P1.pos==0) {
						Hexa.push_back(P1);
						Hexa.push_back(P1.add(0.0f, yStep*get_Coef(P1, P1.add(0.0f, yStep))));
						Hexa.push_back(P2.add(-xStep*get_Coef(P2, P2.add(-xStep, 0.0f)), 0.0f));
						Hexa.push_back(P2);
						Hexa.push_back(P2.add(0.0f, -yStep*get_Coef(P2, P2.add(0.0f, -yStep))));
						Hexa.push_back(P1.add(xStep*get_Coef(P1, P1.add(xStep, 0.0f)), 0.0f));
					}
					if (P1.pos==1) {
						Hexa.push_back(P1);
						Hexa.push_back(P1.add(xStep*get_Coef(P1, P1.add(xStep, 0.0f)), 0.0f));
						Hexa.push_back(P2.add(0.0f, yStep*get_Coef(P2, P2.add(0.0f, yStep))));
						Hexa.push_back(P2);
						Hexa.push_back(P2.add(-xStep*get_Coef(P2, P2.add(-xStep, 0.0f)), 0.0f));
						Hexa.push_back(P1.add(0.0f, -yStep*get_Coef(P1, P1.add(0.0f, -yStep))));
					}

					// On affiche l'hexagone
					glBegin(GL_POLYGON);
						glVertex2fv(Hexa.at(0).getVertex());
						glVertex2fv(Hexa.at(1).getVertex());
						glVertex2fv(Hexa.at(2).getVertex());
						glVertex2fv(Hexa.at(3).getVertex());
						glVertex2fv(Hexa.at(4).getVertex());
						glVertex2fv(Hexa.at(5).getVertex());
					glEnd();
				}
			}
			else {

				// sinon -> Affichage brut : On relie les deux sommets
				glBegin(GL_LINES);
					glVertex2fv(Sq->at(0).getVertex());
					glVertex2fv(Sq->at(1).getVertex());
				glEnd();
			}
		}

		// Si le carré à trois sommets actifs ********************************************************
		else if (Sq->size()==3) {

			// et si on n'est pas en affichage brut
			if (dis_mode > 0) {
				Point P1 = Sq->at(0);
				Point P2 = Sq->at(1);
				Point P3 = Sq->at(2);
				vector<Point> Penta;
				Penta.clear();

				// On calcule les sommets du pentagone à afficher
				if ((P1.pos==0)&&(P2.pos==1)&&(P3.pos==2)) {
					Penta.push_back(P1);
					Penta.push_back(P2);
					Penta.push_back(P3);
					Penta.push_back(P3.add(0.0f, -yStep*get_Coef(P3, P3.add(0.0f, -yStep))));
					Penta.push_back(P1.add(xStep*get_Coef(P1, P1.add(xStep, 0.0f)), 0.0f));
				}
				if ((P1.pos==0)&&(P2.pos==1)&&(P3.pos==3)) {
					Penta.push_back(P1);
					Penta.push_back(P2);
					Penta.push_back(P2.add(xStep*get_Coef(P2, P2.add(xStep, 0.0f)), 0.0f));
					Penta.push_back(P3.add(0.0f, yStep*get_Coef(P3, P3.add(0.0f, yStep))));
					Penta.push_back(P3);
				}
				if ((P1.pos==0)&&(P2.pos==2)&&(P3.pos==3)) {
					Penta.push_back(P1);
					Penta.push_back(P1.add(0.0f, yStep*get_Coef(P1, P1.add(0.0f, yStep))));
					Penta.push_back(P2.add(-xStep*get_Coef(P2, P2.add(-xStep, 0.0f)), 0.0f));
					Penta.push_back(P2);
					Penta.push_back(P3);
				}
				if ((P1.pos==1)&&(P2.pos==2)&&(P3.pos==3)) {
					Penta.push_back(P1);
					Penta.push_back(P2);
					Penta.push_back(P3);
					Penta.push_back(P3.add(-xStep*get_Coef(P3, P3.add(-xStep, 0.0f)), 0.0f));
					Penta.push_back(P1.add(0.0f, -yStep*get_Coef(P1, P1.add(0.0f, -yStep))));
				}

				// On affiche le pentagone
				glBegin(GL_POLYGON);
					glVertex2fv(Penta.at(0).getVertex());
					glVertex2fv(Penta.at(1).getVertex());
					glVertex2fv(Penta.at(2).getVertex());
					glVertex2fv(Penta.at(3).getVertex());
					glVertex2fv(Penta.at(4).getVertex());
				glEnd();
			}
			else {

				// sinon -> Affichage brut : On trace le triangle reliant les trois points
				glBegin(GL_TRIANGLES);
					glVertex2fv(Sq->at(0).getVertex());
					glVertex2fv(Sq->at(1).getVertex());
					glVertex2fv(Sq->at(2).getVertex());
				glEnd();
			}
		}

		// si le carré à ses quatre sommets actifs ******************************************************
		else if (Sq->size()==4) {

			// On trace le carré complet
			glBegin(GL_QUADS);
				glVertex2fv(Sq->at(0).getVertex());
				glVertex2fv(Sq->at(1).getVertex());
				glVertex2fv(Sq->at(2).getVertex());
				glVertex2fv(Sq->at(3).getVertex());
			glEnd();
		}
	}
}

/* Méthode d'affichage des Blobs */
void displayBlobs()
{
	Point P;
	vector<Point> crnt_Square;
	VtoShow.clear();

	// Si l'utilisateur a saisi quelque chose
	if (N > 0) {

		// Pour chaque position de la discrétisation de la fenêtre d'affichage : recensement des carrés à afficher
		for (int i=0; i<maxPt; i++) {
			for (int j=0; j<maxPt; j++) {

				// (Ré-)Initialisation
				crnt_Square.clear();

				/* Examen des 4 sommets dans le sens trigonométrique pour l'affichage OpenGL sans culling */
				P = Point(i*xStep, j*yStep);
				if (is_Inside(P)) {						
					P.setPos(0);
					crnt_Square.push_back(P);						
				}
				P = Point(i*xStep, (j+1)*yStep);
				if (is_Inside(P)) {						
					P.setPos(1);
					crnt_Square.push_back(P);						
				}
				P = Point((i+1)*xStep, (j+1)*yStep);
				if (is_Inside(P)) {
					P.setPos(2);
					crnt_Square.push_back(P);						
				}
				P = Point((i+1)*xStep, j*yStep);
				if (is_Inside(P)) {
					P.setPos(3);
					crnt_Square.push_back(P);						
				}
				// Si le carré considéré a au moins un sommet actif, on l'ajoute à l'affichage
				if (crnt_Square.size()>0) VtoShow.push_back(crnt_Square);
			}
		}
	}
	// On affiche les carrés recensés
	displaySquares();
}


/* Méthode d'affichage GLUT */ 
void main_display(void)
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glShadeModel(GL_SMOOTH);
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);

 	xStep = width/maxPt;		// recalcul du pas de discrétisation horizontal
	yStep = height/maxPt;		// recalcul du pas de discrétisation vertical
	
	glLoadIdentity();
	glPushMatrix();

    	// Affichage des Points saisis par l'utilisateur
		glInitNames();
		glPushName(1);
		for (int i=0; i<N ;i++) {
			glLoadName(i);
			glColor3f(0.0f, 0.96f, 0.0f);
			glPointSize(7.0f);
			glBegin(GL_POINTS);
				glVertex2fv(Pts[i].getVertex());
			glEnd();
		}

		// Affichge des Blobs
		displayBlobs();			

	glPopMatrix();
	
	/* Nettoyage */
	glFlush();
    glutSwapBuffers();
}

/* Méthode GLUT de gestion du clic de la souris */
void Mouse(int button, int state, int x, int y)
{
	GLint viewport[4];
	glutSetCursor(GLUT_CURSOR_CROSSHAIR);
	glGetIntegerv(GL_VIEWPORT, viewport);
	
	// Le bouton gauche sert à saisir de nouveaux points dans la fenêtre
	if ((button == GLUT_LEFT_BUTTON) && (N < NMAX)) {
		droite = false;
		gauche = true;
		glColor3f(0.0f, 0.96f, 0.0f);
		glPointSize(7.0f);
		glInitNames();
		glPushName(1);
		Pts[N].set((float)x, (float)(height-y));
		glLoadName(N);
		glBegin(GL_POINTS);
			glVertex2fv(Pts[N].getVertex());
		glEnd();
		if(state == GLUT_UP) N++;
		glutPostRedisplay();
	}

	// Le bouton droit sert à selectionner le point à déplacer
	if(button == GLUT_RIGHT_BUTTON) {
		gauche = false;
		droite = true;
		if(state == GLUT_DOWN) {
			GLuint *selectBuf = new GLuint[2*NMAX];
			GLuint *ptr;
			GLint hits;
			glSelectBuffer(NMAX, selectBuf);
			glRenderMode(GL_SELECT);
			glPushMatrix();
				glMatrixMode(GL_PROJECTION);
				glLoadIdentity();
				gluPickMatrix((float)x, (float)(height-y), 14.0f, 14.0f, viewport);
				glOrtho(0.0f, (float)width, 0.0f, (float)height, -1.0f, 1.0f);
				glColor3f(0.0f, 0.96f, 0.0f);
				glPointSize(7.0f);
				glInitNames();
				glPushName(1);
				for (int i=0; i<N; i++) {
					glLoadName(i);
					glBegin(GL_POINTS);
						glVertex2fv(Pts[i].getVertex());
					glEnd();			
				}
			glPopMatrix();
			glFlush();
			hits = glRenderMode(GL_RENDER);
			if(hits) {
				ptr = (GLuint *)selectBuf;
				ptr += 3;
				selected = *ptr;
			}
		}
		if(state == GLUT_UP) selected = -1;
		main_reshape(viewport[2], viewport[3]);
		glutPostRedisplay();
	}
}

/* Méthode GLUT de gestion du mouvement de la souris */
void Motion(int xi, int yi)
{
	float x = (float)xi;
	float y = (float)yi;
	if ((droite) && (selected!=-1)) {

		// On contraint le curseur dans la fenêtre GLUT
		if (x<10.0f) x=10.0f;
		if (x>width-10.0f) x=width-10.0f;
		if (y<10.0f) y=10.0f;
		if (y>height-10.0f) y=height-10.0f;

		// On modifie les coordonnées du point selectionné
		int i= selected;
		Pts[i].set(x, (float)height-y);
		glutPostRedisplay();
	}
}

/* Méthode GLUT de gestion des saisies clavier */
void main_keyboard(unsigned char c, int x, int y)
{
	switch (c) {
		
		case 'S':	if (seuil<1.0f) seuil += 0.1f;	break;		// On augmente la valeur du seuil
		case 's':	if (seuil>0.1f) seuil -= 0.1f;	break;		// On diminue la valeur du seuil

		case 'm':	mix_mode = 0;	break;		// mélange des potentiels
		case 'u':	mix_mode = 1;	break;		// intersection des potentiels
		case 'i':	mix_mode = -1;	break;		// union des potentiels

		case 'M':	pot_mode = 1;	break;		// potentiel de Murakami
		case 'N':	pot_mode = 2;	break;		// potentiel de Nishimura
		case 'W':	pot_mode = 3;	break;		// potentiel de Wyvill
		case 'B':	pot_mode = 0;	break;		// potentiel de Blinn (a=1 et b=1)

		case 'z':	dis_mode = 0;	break;		// affichage brut
		case 'd':	dis_mode = 1;	break;		// affichage par les milieux
		case 'p':	dis_mode = 2;	break;		// affichage proportionnel

		default :	break;
	}
	glutPostRedisplay();
}

/* Lanceur de l'application avec initialisation GLUT */
int main (int argc, char** argv)
{
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_DEPTH | GLUT_DOUBLE);
    glutInitWindowPosition(0, 0);
    glutInitWindowSize(width, height);
    glutCreateWindow("Select");
    glEnable(GL_DEPTH_TEST);

    glutDisplayFunc(main_display);
	glutKeyboardFunc(main_keyboard);
  	glutMouseFunc(Mouse);
	glutMotionFunc(Motion);	
	glutReshapeFunc(main_reshape);

    glutMainLoop();
    return 0;
}