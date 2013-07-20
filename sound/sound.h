/* Copyright 2013 by trilliondimension@163.com */

#ifndef SOUND_H
#define SOUND_H

#include <SDL/SDL.h>
#include <SDL/SDL_mixer.h>

#include <stdio.h>

class SoundWav
{
private:
	Mix_Chunk *wav;
	int channel;
	int volume;
public:
	SoundWav() : wav(NULL), channel(-1),  volume(100){}
	SoundWav(char *filename);

	~SoundWav();

	void load(char *file);


	// loops  = -1 is infinit
	void play(int loop);


	int isPlaying();


	void stop();


	int setVolume(int level);


	void pause();

	int isPaused();

	void resume();

};

class SoundMusic
{
private:
	Mix_Music *music;
	int volume;
	
public:
	SoundMusic();

	SoundMusic(char *file);
	
	~SoundMusic();

	void load(char *file);


	void play(int loop);


	void stop();

	int setVolume(int level);

	int isPlaying();

	void resume();
	
	void pause();
	
	int isPaused();


	void rewind();

};


#endif
