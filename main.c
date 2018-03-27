#include "portaudio.h"


typedef struct
{
	float left_phase;
	float right_phase;
}
paTestData;

static int patestCallback(	const void *inputBuffer, void *outputBuffer,
			   	unsigned long framesPerBuffer, 
				const PaStreamCallbackTimeInfo* timeInfo, 
				PaStreamCallbackFlags statusFlags, 
				void *userData)
{
	paTestData *data = (paTestData*) userData;
	float *out = (float*) outputBuffer;
	unsigned int i;
	(void) inputBuffer;

	for( i = 0 ; i < framesPerBuffer ; i++)
	{
		*out++ = data->left_phase;
		*out++ = data->right_phase;
		data->left_phase += 0.01f;
		if(data->left_phase >= 1.0f) data->left_phase -= 2.0f;
		data->right_phase += 0.03f;
		if( data->right_phase >= 1.0f) data->right_phase -= 2.0f;
	}
	return 0;
}
