#include <GL/glut.h>
#include <unistd.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include "sound.h"

void rect(void)
{
	glClear(GL_COLOR_BUFFER_BIT);
	glColor3f(1.0,1.0,0.0);
	glRectf(-0.25 * 100,0.25 * 100,0.25*100,-0.25*100);
	
	glFlush();
}

void SetupRC()
{
	glClearColor(0.0,0.0,1.0,1.0);
}

void ChangeSize(int width, int height)
{
	if(height == 0)
		height = 1;
	glViewport(0,0,width,height);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	float aspecRatio = (float)width/height;
	if(width <= height)
		glOrtho(-100.0, 100.0, -100/aspecRatio,100.0/aspecRatio,1.0,-1.0);
	else
		glOrtho(-100.0 * aspecRatio, 100.0 * aspecRatio, -100.0, 100.0, 1.0,-1.0);
	
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

int main(int argc, char *argv[])
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	glutCreateWindow("rect");

	glutDisplayFunc(rect);
	glutReshapeFunc(ChangeSize);
	SetupRC();
	
	SDL_Init(SDL_INIT_AUDIO);	

	Mix_Init(MIX_INIT_FLAC | MIX_INIT_MP3 | MIX_INIT_MOD | MIX_INIT_OGG);

	Mix_OpenAudio(44100, AUDIO_S16SYS, 2, 1024);
	
	SoundWav *s = new SoundWav;
	
	s->load("goggles.wav");
	s->play(0);
	int i = 100;
	int flag = 1;
//	while(1)
//	{
//		s->setVolume(i);
//		usleep(80000);
//			i--;
//	}

//	while(s->isPlaying());
//	s->load("thunder1.wav");
//	s->setVolume(10);
//	s->play(3, 0);

	SoundMusic music("youRock.mp3");
	music.play(0);
	music.setVolume(80);
	music.pause();
	sleep(15);
	music.resume();



	glutMainLoop();
	Mix_CloseAudio();
	
	Mix_Quit();
	return 0 ;
}
