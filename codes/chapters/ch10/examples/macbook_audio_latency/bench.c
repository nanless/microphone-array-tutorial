#include <portaudio.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

enum {RATE=48000, FIRST=24000, PERIOD=28800, PROBE=960};
typedef struct {unsigned long count, flags; size_t start; double adc, now, dac;} Event;
typedef struct {
    size_t pos, total, events, capacity;
    int probes, muted, delay, alternating, gate_begin, gate_length;
    float *input, *output;
    float probe[PROBE];
    Event *log;
} Bench;

static int callback(const void *in_, void *out_, unsigned long count,
    const PaStreamCallbackTimeInfo *times, PaStreamCallbackFlags flags, void *data) {
    Bench *b=data;
    const float *in=in_;
    float *out=out_;
    if (b->events>=b->capacity) return paAbort;
    b->log[b->events++]=(Event){count,flags,b->pos,times->inputBufferAdcTime,times->currentTime,times->outputBufferDacTime};
    for(unsigned long j=0;j<count;j++) {
        size_t i=b->pos+j;
        float y=0;
        if(i<b->total) {
            b->input[i]=in?in[j]:0;
            if(i>=FIRST) {
                size_t k=(i-FIRST)%PERIOD;
                size_t trial=(i-FIRST)/PERIOD;
                int delay=b->delay+(b->alternating && trial%2 ? 1920 : 0);
                if(trial<(size_t)b->probes) {
                    if(k<PROBE && !b->muted) y=b->probe[k];
                    // This finite gate prevents feedback. Only the first arrival
                    // is forwarded, exactly `delay` samples later.
                    if(delay && k>=(size_t)(b->gate_begin+delay)
                       && k<(size_t)(b->gate_begin+delay+b->gate_length))
                        y+=fmaxf(-0.14f,fminf(0.14f,0.28f*b->input[i-delay]));
                }
            }
            b->output[i]=y;
        }
        out[2*j]=out[2*j+1]=y;
    }
    b->pos+=count;
    return b->pos>=b->total?paComplete:paContinue;
}

static void check(PaError err) {
    if(err<0) {fprintf(stderr,"%s (%d)\n",Pa_GetErrorText(err),err); exit(2);}
}

static FILE *create(const char *base,const char *suffix,const char *mode) {
    char file[4096];
    if(snprintf(file,sizeof file,"%s.%s",base,suffix)>=(int)sizeof file) exit(2);
    FILE *f=fopen(file,mode);
    if(!f) {perror(file); exit(2);}
    return f;
}

int main(int argc,char **argv) {
    if(argc!=8) {
        fprintf(stderr,"usage: bench prefix frames probes mute delay_samples gate_begin gate_length\n");
        return 2;
    }
    int frames=atoi(argv[2]);
    Bench b={.probes=atoi(argv[3]),.muted=atoi(argv[4]),.delay=atoi(argv[5]),
             .gate_begin=atoi(argv[6]),.gate_length=atoi(argv[7])};
    if(b.delay<0 && b.delay>=-RATE) {b.alternating=1; b.delay=-b.delay;}
    if(frames<32||frames>1024||b.probes<1||b.probes>64||b.delay<0||
       b.gate_begin<0||b.gate_length<0||b.gate_begin+b.gate_length+b.delay+(b.alternating ? 1920 : 0)>=PERIOD||
       (b.delay && b.delay<=b.gate_length)) return 2;
    b.total=FIRST+PERIOD*(size_t)b.probes;
    b.capacity=b.total/16+8;
    b.input=calloc(b.total,sizeof(float)); b.output=calloc(b.total,sizeof(float));
    b.log=calloc(b.capacity,sizeof(Event));
    if(!b.input||!b.output||!b.log) return 2;
    for(int i=0;i<PROBE;i++) {
        double t=i/(double)RATE, envelope=sin(M_PI*i/(PROBE-1));
        b.probe[i]=(float)(0.14*envelope*envelope*sin(2*M_PI*(900*t+(4500-900)*t*t/(2.0*PROBE/RATE))));
    }
    check(Pa_Initialize());
    PaDeviceIndex input=Pa_GetDefaultInputDevice(), output=Pa_GetDefaultOutputDevice();
    if(input==paNoDevice||output==paNoDevice) return 2;
    const PaDeviceInfo *idi=Pa_GetDeviceInfo(input), *odi=Pa_GetDeviceInfo(output);
    PaStreamParameters ip={input,1,paFloat32,idi->defaultLowInputLatency,NULL};
    PaStreamParameters op={output,2,paFloat32,odi->defaultLowOutputLatency,NULL};
    PaStream *stream;
    check(Pa_OpenStream(&stream,&ip,&op,RATE,frames,paNoFlag,callback,&b));
    const PaStreamInfo info=*Pa_GetStreamInfo(stream);
    printf("input=%s\noutput=%s\nhost=%s\nPortAudio=%s\nrate=%.0f\ninput_latency_s=%.9f\noutput_latency_s=%.9f\n",
        idi->name,odi->name,Pa_GetHostApiInfo(idi->hostApi)->name,Pa_GetVersionText(),info.sampleRate,info.inputLatency,info.outputLatency);
    fflush(stdout);
    check(Pa_StartStream(stream));
    int active;
    while((active=Pa_IsStreamActive(stream))==1) Pa_Sleep(25);
    check(active); check(Pa_CloseStream(stream)); check(Pa_Terminate());
    if(b.pos<b.total) return 3;
    FILE *log=create(argv[1],"csv","w");
    fprintf(log,"start,frames,flags,adc_s,callback_s,dac_s\n");
    for(size_t i=0;i<b.events;i++){
        Event e=b.log[i];
        fprintf(log,"%zu,%lu,%lu,%.12f,%.12f,%.12f\n",e.start,e.count,e.flags,e.adc,e.now,e.dac);
    }
    fclose(log);
    FILE *p=create(argv[1],"probe.f32","wb"); fwrite(b.probe,sizeof(float),PROBE,p); fclose(p);
    p=create(argv[1],"input.f32","wb"); fwrite(b.input,sizeof(float),b.total,p); fclose(p);
    p=create(argv[1],"output.f32","wb"); fwrite(b.output,sizeof(float),b.total,p); fclose(p);
    printf("callbacks=%zu\n",b.events);
    free(b.input); free(b.output); free(b.log);
}
