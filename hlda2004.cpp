#include "mex.h"
#include "matrix.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "lightspeed\util.h"
#include <time.h>


#define PARAM_MAXITER 10000
#define PARAM_INITITER 50



// Global variables
int g_ndocs=0;
double* g_eta;
double* g_alpha;
int g_depth;
double* g_gamma;
int g_nterms;


clock_t time1,time2,time3,time4;
double diff1,diff2,diff3,diff4;

typedef struct _topic{
	int* words; // Frequency matrix of topic
	struct _topic* next; // Linked list 
	struct _topic* prev; // Double Linked List
	struct _topic* parent; // Parent node
	struct _topic* copy; // Used to facilitate copy of topics
	int ndocs;
	int level;
	int nwords;
	int id;
	double loglikelihood; // Gamma ratio
} Topic;


typedef struct _topiclist{
	Topic* head;
	Topic* tail;
	int count;
} TopicList;

typedef struct _doc{
	int* words;
	int* levels;
	int* Bow;  // Bag of Words
	int* indexset; // Bow is sparse
	int ndiffwords;
	int length;
	Topic** path;
	int* levelsum;
} Doc;



void addtoTail(TopicList* tl,Topic* t)
{
		tl->tail->next = t;
		t->prev = tl->tail;
		tl->tail = t;
		tl->count++;
}

void removeTopic(TopicList* tl,Topic* t)
{

	if (t == tl->head) 
	{
		if (tl->head != tl->tail)
		{
			tl->head = t->next;
			tl->head->prev = NULL;
		}
		else
		{
			tl->head = tl->tail = NULL;
		}
		
	}
	else if(t == tl->tail)
	{
		if (tl->head != tl->tail)
		{
			tl->tail = t->prev;
			tl->tail->next = NULL;
		}
		else
		{
			tl->head = tl->tail = NULL;
		}
	}
	else
	{
		t->prev->next = t->next;
		t->next->prev = t->prev;
	}
	mxFree(t->words);
	// mxFree(t); // Batch allocation could not be free'd indivudually
	tl->count--;
}

Topic* emptyTopic()
{
	Topic* t = (Topic*) mxCalloc(1,sizeof(Topic));
	t->level = 0;
	t->ndocs = 0;
	t->next = NULL;
	t->prev = NULL;
	t->parent = NULL;
	t->copy = NULL;
	t->words = (int*) mxCalloc(g_nterms,sizeof(int));
	t->nwords = 0;
	return t;
}




// Dr. Al Hasan's Algorithm lecture
void inPlacePerm(int* data,int size)
{
	int i,temp,r;
	for (i=0;i<size;i++)
		{
			r = i + ( rand() % (size-i) );
			temp = data[i];
			data[i] = data[r];
			data[r] = temp;
		}
}



// C uniform random variable approximation 
double urand()
{
	return rand() / ((double) (RAND_MAX+1.0));
}

Doc* readfile(char* filename)
{
	FILE* file;
	int doccount=0;
	Doc* docs;

	int n,word,count;
	int w,i,j,k;
	file = fopen(filename,"r");
	if (file==0)
	{
		printf("Unable to open file\n");
	}
	else
	{
		// Format sample : 9 (number of words)  9(word):6(counts) 20:4
		// Expands documents into plain format
		// ! hz ! Later I can add support to plain format , even folder scan for documents
		while (!feof(file)&& (fscanf(file,"%d",&n)>0))
		{
			doccount++;
			for(w=0;w<n;w++)
			{
				fscanf(file,"%d:%d",&word,&count);
			}
		}

		fseek(file, 0L, SEEK_SET);
		docs = (Doc*) mxCalloc(doccount,sizeof(Doc));

		// Get length of files
		for (i=0;i<doccount;i++)
		{
			docs[i].length = 0;
			fscanf(file,"%d",&n); // How many words exist in doc
			for(w=0;w<n;w++)
			{
				fscanf(file,"%d:%d",&word,&count);
				docs[i].length+=count;
			}
		}

		fseek(file, 0L, SEEK_SET);
		// Get actual words , expand it
		for (i=0;i<doccount;i++)
		{
			k=0;
			fscanf(file,"%d",&n); // How many words exist in doc

		
			// Initial values of docs
			docs[i].levels = (int*) mxCalloc(docs[i].length,sizeof(int));
			docs[i].words = (int*) mxCalloc(docs[i].length,sizeof(int));
			docs[i].levelsum = (int*) mxCalloc(g_depth,sizeof(int));
			docs[i].path = (Topic**) mxCalloc(g_depth,sizeof(Topic*));
			docs[i].Bow =  (int*) mxCalloc(g_depth*g_nterms,sizeof(int));
			docs[i].indexset = (int*) mxCalloc(n,sizeof(int));
			docs[i].ndiffwords = n;

			for(w=0;w<n;w++)
			{
				fscanf(file,"%d:%d",&word,&count);
				docs[i].indexset[w] = word;
				for(j=0;j<count;j++,k++)
				{
					docs[i].levels[k] = 0;
					docs[i].words[k]  = word;
				}

			}

			inPlacePerm(docs[i].words,docs[i].length);
		}

		fclose(file);
	}
	g_ndocs = doccount;
	return docs;
}


TopicList* copyTopics(TopicList* topics)
{
	TopicList* copy;
	Topic* temp;
	int i,j;

	copy = (TopicList*) mxCalloc(1,sizeof(TopicList));
	copy->count = topics->count;

	copy->head = (Topic*) mxCalloc(copy->count,sizeof(Topic));
	copy->tail = copy->head + copy->count - 1; // Bulk allocation

	for(i=0,temp=topics->head;i<topics->count;i++,temp=temp->next)
	{
		temp->copy = copy->head + i;
		copy->head[i].level = temp->level;
		copy->head[i].loglikelihood = temp->loglikelihood;
		copy->head[i].ndocs = temp->ndocs;
		copy->head[i].nwords = temp->nwords;
		copy->head[i].next = copy->head + i +1;
		copy->head[i].prev = copy->head + i -1;
		
		if (i!=0) //Not Root
			copy->head[i].parent = temp->parent->copy; //Parent is always added first

		copy->head[i].words = (int*) mxCalloc(g_nterms,sizeof(int));
		for(j=0;j<g_nterms;j++)
			copy->head[i].words[j] = temp->words[j];
	}

	return copy;
}

void  freeDocs(Doc* docs)
{
	int i;
	for (i=0;i<g_ndocs;i++)
	{
		mxFree(docs[i].levels);
		mxFree(docs[i].words);
		mxFree(docs[i].levelsum);
		mxFree(docs[i].path);
		mxFree(docs[i].Bow);
	}
	mxFree(docs);
}

Doc* copyDocs(Doc* docs)
{
	int i,j,l;
	Doc* copy;

	copy = (Doc*) mxCalloc(g_ndocs,sizeof(Doc));

	for(i=0;i<g_ndocs;i++)
	{
		copy[i].length = docs[i].length;
		copy[i].ndiffwords = docs[i].ndiffwords;
		copy[i].Bow = (int*) mxCalloc(g_nterms*g_depth,sizeof(int));
		copy[i].levelsum = (int*) mxCalloc(g_depth,sizeof(int));
		copy[i].path = (Topic**) mxCalloc(g_depth,sizeof(Topic*));
		copy[i].levels = (int*) mxCalloc(copy[i].length,sizeof(int));
		copy[i].words = (int*) mxCalloc(copy[i].length,sizeof(int));
		copy[i].indexset = (int*) mxCalloc(copy[i].ndiffwords,sizeof(int));

		for(j=0;j<copy[i].length;j++)
		{
			copy[i].levels[j] = docs[i].levels[j];
			copy[i].words[j] = docs[i].words[j];
		}

		for(j=0;j<g_depth;j++)
		{
			copy[i].path[j] = docs[i].path[j]->copy; // Topics need to be copied first
			copy[i].levelsum[j] = docs[i].levelsum[j];
		}

		for(j=0;j<copy[i].ndiffwords;j++)
			copy[i].indexset[j] = docs[i].indexset[j];

		for(j=0;j<g_nterms;j++)
			for(l=0;l<g_depth;l++)
				copy[i].Bow[l+j*g_depth] = docs[i].Bow[l+j*g_depth];
	}
	return copy;
}

void addWords(Doc* doc,Topic* topic)
{
	int i,j,level;
	level = topic->level;
	for (j=0;j<doc->ndiffwords;j++)
	{
		i = doc->indexset[j];
		topic->words[i] += doc->Bow[ level  + g_depth*i];
	}
	topic->nwords += doc->levelsum[level]; 
}

void removeWords(Doc* doc,Topic* topic)
{
	int i,j,level;
	level = topic->level;
	for (j=0;j<doc->ndiffwords;j++)
	{
		i = doc->indexset[j];
		topic->words[i] -= doc->Bow[ level  + g_depth*i];
	}
	topic->nwords -= doc->levelsum[level]; 
}


void addToPath(Doc* doc,Topic* topic,TopicList* tl)
{
	int i,j;
	Topic *temp,*prev;


	for (i=topic->level,temp=topic;i>=0;i--)
	{
		doc->path[temp->level] = temp;
		temp->ndocs++;
		addWords(doc,temp);
		temp = temp->parent;
	}

	// Create new topics  if selected topic is internal node
	for (i=topic->level+1,prev=topic;i<g_depth;i++)
	{
		temp = emptyTopic();
		temp->level  = i;
		temp->ndocs++;
		temp->parent = prev;
		addtoTail(tl,temp);
		doc->path[temp->level] = temp;
		addWords(doc,temp);
		prev = temp;
	}

}

void removeFromPath(Doc* doc,TopicList* topics)
{
	int i,j;
	for (i=g_depth-1;i>=0;i--) // From leaf to parent
	{
		doc->path[i]->ndocs--;
		removeWords(doc,doc->path[i]);
		//Remove unnecessary nodes
		if (doc->path[i]->ndocs == 0 && i!=0) 
		{
				removeTopic(topics,doc->path[i]);
		}
		doc->path[i] = NULL;
	}
	
	


}


void sampleLevels(Doc* doc,int init)
{
	int i,j,word,level;
	double* topicpresence = (double*) mxCalloc(g_depth,sizeof(double));
	double* topiclikelihood = (double*) mxCalloc(g_depth,sizeof(double));
	double* prob = (double*) mxCalloc(g_depth,sizeof(double));
	double  denominator,arand;

	for(i=0;i<doc->length;i++)
	{
		word = doc->words[i];
		level = doc->levels[i];

		if (init != 1)
		{
			// Remove from topic
			doc->path[level]->words[word]--;
			doc->path[level]->nwords--;
			doc->levelsum[level]--;
			doc->Bow[level + g_depth*word ]--;
		}

		// Gibbs probabilities
		denominator=0;
		for (j=0;j<g_depth;j++)
		{
			topicpresence[j] = doc->levelsum[j] + g_alpha[j];
			topiclikelihood[j] = ( doc->path[j]->words[word] + g_eta[j] ) / ( doc->path[j]->nwords + g_nterms*g_eta[j]) ;
			prob[j] = topicpresence[j] * topiclikelihood[j];
			denominator+=prob[j];
		}

		//Sample
		arand = urand() * denominator;
		for(j=0;j<g_depth;j++)
		{
			if 	(arand <= prob[j])
			{
				level = j;
				break;
			}
			else
				arand = arand - prob[j];
		}

		// Update
		doc->levelsum[level]++;
		doc->path[level]->words[word]++;
		doc->path[level]->nwords++;
		doc->levels[i] = level;
		doc->Bow[level + g_depth*word ]++;


	}


	mxFree(topicpresence);
	mxFree(topiclikelihood);
	mxFree(prob);

}

// If new empty topics created after internal node what is likelihood of non created topics
double* emptyTopicScore(Doc* doc)
{
	int w,i,l;
	double* logsum,*cumlogsum;

	logsum = (double*) mxCalloc(g_depth,sizeof(double));
	cumlogsum = (double*) mxCalloc(g_depth,sizeof(double));


	// Same as gammaRatio with empty topic
	for(l=0;l<g_depth;l++)
	{
		logsum[l] = 0;
		for (i=0;i<doc->ndiffwords;i++)
		{
			w = doc->indexset[i];
			if (doc->Bow[l+g_depth*w]==0) continue; // numerator and denominator is same
			logsum[l] += gammaln(doc->Bow[l+g_depth*w] + g_eta[l]) - gammaln(g_eta[l]);
		}
		logsum[l] += gammaln(g_nterms*g_eta[l]);
		logsum[l] -= gammaln(g_nterms*g_eta[l]+ doc->levelsum[l]);
	}
	

	for(l=0;l<g_depth;l++)
	{
		cumlogsum[l] = 0;
		for(i=g_depth-1;i>l;i--) // How many empty nodes required to fill the path
			cumlogsum[l] += logsum[i];
	}
	logsum;
	mxFree(logsum);
	return cumlogsum;
}

double gammaRatio(Doc* doc,Topic* t)
{
	int w,j,l;
	double logsum;

	l = t->level;
	logsum= 0;
	for (j=0;j<doc->ndiffwords;j++)
	{
		w = doc->indexset[j]; // Do not iterate the words that has zero count , basically because ratio is 1 
		if (doc->Bow[l+g_depth*w]==0) continue; // Also for this level
		

		logsum += gammaln(doc->Bow[l+g_depth*w] + g_eta[l] + t->words[w]);
		logsum -= gammaln(g_eta[l] + t->words[w]);
	}

	logsum += gammaln(g_nterms*g_eta[l]+t->nwords);
	logsum -= gammaln(g_nterms*g_eta[l]+ doc->levelsum[l]+t->nwords);

	return logsum;
}



void samplePath(Doc* doc,TopicList* topics,int init)
{
	double *logempty,*ncrp,*likelihood,*gammaratio,*probs;  // All in logarithm
	double maxprob,denominator,arand,pval;
	Topic *topic,*temp;
	int i,j;

	likelihood = (double*) mxCalloc(topics->count,sizeof(double));
	ncrp  = (double*) mxCalloc(topics->count,sizeof(double));
	probs = (double*) mxCalloc(topics->count,sizeof(double));


	
	if (!init)
	{
		removeFromPath(doc,topics);
	}

	logempty = emptyTopicScore(doc); // Cumulative
	

	// Calculate likelihood value and do not iterate many times on parents
	
	for(i=0,topic = topics->head;i<topics->count;i++)
	{
		topic->loglikelihood = gammaRatio(doc,topic);
		topic = topic->next;
	}
	

	for(i=0,topic = topics->head;i<topics->count;i++)
	{

		
		

		// nCRP
		ncrp[i] = 0;
		temp = topic;
		for(j=topic->level;j>0;j--)
		{
			ncrp[i] += log( temp->ndocs / ((temp->parent)->ndocs + g_gamma[j-1] )); // Parent's level gamma
			temp = temp->parent;
		}

		if (topic->level<g_depth-1) // Internal node
			ncrp[i] += log( g_gamma[topic->level] / (topic->ndocs + g_gamma[topic->level] ) );

		
		

		//Likelihood
		likelihood[i] =0;
		temp = topic;
		for(j=topic->level;j>=0;j--)
		{
			likelihood[i] += temp->loglikelihood; 
			temp = temp->parent;
		}

		likelihood[i] += logempty[topic->level];
		probs[i] = likelihood[i] + ncrp[i];
		topic = topic->next;

		
	}
	


	// SAMPLE TOPIC
	denominator=0;
	maxprob=probs[0];
	arand = urand();

	//Avoid numerical underflow
	for(i=0;i<topics->count;i++){ 
		if (maxprob<probs[i])
			maxprob = probs[i];
	}

	for(i=0;i<topics->count;i++){ 
		probs[i] -= maxprob;
		probs[i] = exp(probs[i]); // Convert to actual prob
		denominator += probs[i];
	}
	
	//Select Topic
	for(i=0,topic = topics->head;i<topics->count;i++)
	{ 
		pval = (probs[i]/denominator);
		if (pval>arand)
		{
			temp = topic; // Selected topic
			break;
		}
		else
		{
			arand = arand - pval;
		}
		topic=topic->next;
	}
	


	// Now add to tree
	addToPath(doc,temp,topics);

	mxFree(ncrp);
	mxFree(likelihood);
	mxFree(logempty);
	mxFree(probs);
}


double etaScore(TopicList* topics)
{
	
	int i,j;
	double logsum;
	Topic* t;
	
		
	logsum= 0;
	for(t=topics->head,j=0;j<topics->count;j++,t=t->next)
	{
		for (i=0;i<g_nterms;i++)
		{
			if (t->words[i]==0) continue;
			logsum += gammaln(g_eta[t->level] + t->words[i]);
			logsum -= gammaln(g_eta[t->level]);

		}

		logsum += gammaln(g_nterms*g_eta[t->level]);
		logsum -= gammaln(g_nterms*g_eta[t->level]+t->nwords);
	}
	return logsum;
}


double alphaScore(Doc* docs)
{
	int i,j;
	double logsum;
	double alphasum,sumgammalnalpha;

	alphasum = 0;
	sumgammalnalpha = 0;
	for(j=0;j<g_depth;j++)
	{
		alphasum += g_alpha[j];
		sumgammalnalpha += gammaln(g_alpha[j]);
	}


	logsum = 0;
	for(i=0;i<g_ndocs;i++)
	{
		logsum += gammaln(alphasum) - sumgammalnalpha - gammaln(docs[i].length+alphasum);
		for(j=0;j<g_depth;j++){
			logsum += gammaln(docs[i].levelsum[j]+g_alpha[j]);
		}
	}
	return logsum;
}



// !! hz !! Check this function 3 times
// Blei has different CRP score (gammaScore)
double crpScore(TopicList* topics)
{
	int i,j;
	double logsum;
	Topic* t;
	FILE* fp;
	
		
	logsum= gammaln(g_gamma[0])  - gammaln(g_gamma[0] + topics->head->ndocs); // Root level

	for(t=topics->head->next,j=1;j<topics->count;j++,t=t->next)
	{
		if (t->level < g_depth-1) // Non leaf nodes
		{
			logsum +=  gammaln(g_gamma[t->level])  - gammaln(g_gamma[t->level] + t->ndocs); // Denominator of crp in this level
		}
		logsum += log(g_gamma[t->parent->level]) + gammaln(t->ndocs); // crp numerator in parents level
	}
	
}


TopicList* reset(Doc* docs)
{
		TopicList* topics;
		int i,j,w,l;
		// Initial Tree
		topics = (TopicList*) mxCalloc(1,sizeof(TopicList));
		topics->head = topics->tail = emptyTopic();
		topics->count =1;

			

		for (i=0;i<g_ndocs;i++)
		{
			for(l=0;l<g_depth;l++)
			{
				docs[i].path[l] = 0;
				docs[i].levelsum[l]=0;
				for(j=0;j<docs[i].ndiffwords;j++)
				{
					w = docs[i].indexset[j];
					docs[i].Bow[l+w*g_depth] = 0; // Others are not used
				}
			}
			for(j=0;j<docs[i].length;j++)
			{
				docs[i].levels[j]  = 0;
			}
			inPlacePerm(docs[i].words,docs[i].length);
		}

		return topics;
}

void halidInitial(Doc* docs,TopicList* topics)
{
	int i,j;

	for (i=0;i<g_ndocs;i++)
		for(j=0;j<docs[i].length;j++)
		{
			docs[i].levels[j] = (rand() % g_depth); // Level starts from 0
			docs[i].Bow[ docs[i].levels[j]+g_depth*docs[i].words[j]]++;
			docs[i].levelsum[docs[i].levels[j]]++;
		}
	for (i=0;i<g_ndocs;i++)
	{	
		samplePath(&docs[i],topics,1);
		sampleLevels(&docs[i],0);
	}

}


void halidMultiInitial(Doc* docs,TopicList* topics)
{
	int i,j;

	for (i=0;i<g_ndocs;i++)
		for(j=0;j<docs[i].length;j++)
		{
			docs[i].levels[j] = (rand() % g_depth); // Level starts from 0
			docs[i].Bow[ docs[i].levels[j]+g_depth*docs[i].words[j]]++;
			docs[i].levelsum[docs[i].levels[j]]++;
		}
	for (i=0;i<g_ndocs;i++)
	{	
		samplePath(&docs[i],topics,1);
		sampleLevels(&docs[i],0);
	}
	for (j=0;j<200;j++)
	{
		for (i=0;i<g_ndocs;i++)
		{	
			samplePath(&docs[i],topics,0);
			sampleLevels(&docs[i],0);
		}
	}

}

void DundarSinglePathInitial(Doc* docs,TopicList* topics)
{
	int i,j;

	// Random sample levels
	for (i=0;i<g_ndocs;i++)
		for(j=0;j<docs[i].length;j++)
		{
			docs[i].levels[j] = (rand() % g_depth); // Level starts from 0
			docs[i].Bow[ docs[i].levels[j]+g_depth*docs[i].words[j]]++;
			docs[i].levelsum[docs[i].levels[j]]++;
		}

	addToPath(&docs[0],topics->head,topics); // Add to root
	for (i=1;i<g_ndocs;i++)
	{	
		addToPath(&docs[i],topics->tail,topics); // Same path
		sampleLevels(&docs[i],0);
	}

}

double* normalize(double* vect)
{
	int i;
	double sum=0;
	double* normalized = (double*) mxCalloc(g_depth,sizeof(double));
	for (i=0;i<g_depth;i++)
		sum += vect[i];
	for (i=0;i<g_depth;i++)
		normalized[i] = vect[i]/sum;

	return normalized;
}

void halidAlphaInitial(Doc* docs,TopicList* topics)
{
	int i,j,k;
	double r;
	double* normalized_alpha = normalize(g_alpha);
	for (i=0;i<g_ndocs;i++)
		for(j=0;j<docs[i].length;j++)
		{
			r = urand();
			for (k=0;k<g_depth;k++)
				if ((r=r-normalized_alpha[k])<0)
					break;
			docs[i].levels[j] = k;
			docs[i].Bow[ docs[i].levels[j]+g_depth*docs[i].words[j]]++;
			docs[i].levelsum[docs[i].levels[j]]++;
		}
	for (i=0;i<g_ndocs;i++)
	{	
		samplePath(&docs[i],topics,1);
		sampleLevels(&docs[i],0);
	}
	mxFree (normalized_alpha);
}


void BleiInitial(Doc* docs,TopicList* topics)
{  // !! hz !! I can have some mistake , it is not working sometimes
	
	int i,j;
	
	for (i=0;i<g_ndocs;i++)
	{
		addToPath(&docs[i],topics->head,topics);
		sampleLevels(&docs[i],1);
		if (i>0) samplePath(&docs[i],topics,0);
		sampleLevels(&docs[i],0);
	
	}
}


mxArray* getParentMatrix(TopicList* topics)
{
	mxArray* m;
	double* out;
	Topic* current,*parent;
	int parentindex,ndocs,i,j;

	m=mxCreateDoubleMatrix(2,topics->count,mxREAL);
	out = mxGetPr(m);

	current = topics->head;
		
	for(i=0;i<topics->count;i++,current=current->next)
	{
			
		parent=topics->head;
		current->id = i+1; //Matlab index

		for (j=0;j<topics->count;j++,parent=parent->next)
			if (current->parent == parent)
			{
				parentindex = j+1; // Matlab index
				ndocs = current->ndocs;
			}
		out[2*i] = parentindex;
		out[2*i+1] = ndocs;
	}
	out[0] = 0; // Root is 0
	out[1] = topics->head->ndocs;
	return m;
}

// One need to call getParentMatrix first
mxArray* getPaths(Doc* docs)
{
	int i,j;
	mxArray* m;
	double* mm ;

	m = mxCreateDoubleMatrix(g_ndocs,g_depth,mxREAL);
	mm = mxGetPr(m);
	for (i=0;i<g_ndocs;i++)
	{
		for(j=0;j<g_depth;j++)
		{
			mm[i+j*g_ndocs]= docs[i].path[j]->id;
		}
	}
	return m;

}

mxArray* getTopics(TopicList* topics)
{
	Topic* current;
	int i,j;
	mxArray* topicmat;
	double* values;

	topicmat = mxCreateDoubleMatrix(topics->count,g_nterms,mxREAL);
	values = mxGetPr(topicmat);

	for (i=0,current=topics->head;i<topics->count;i++,current=current->next)
	{
		for(j=1;j<g_nterms;j++)
			values[i+(j-1)*topics->count] = current->words[(j-1)]; // Back to Matlab index
	}
	return topicmat;
}



void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
                 const mxArray *prhs[])

{


	char* filename;
	FILE* file;
	Doc* docs,*bestdocs;
	TopicList* topics,*besttopics;
	int i,j,iter;
	Topic* root;
	double score,bestscore,lasttime,*scorearray;
	Topic* temp;


		if (nrhs>5)
		{
			filename = (char*) mxCalloc((mxGetN(prhs[0])+1),sizeof(char));

			g_depth = mxGetScalar(prhs[1]);
			g_nterms = mxGetScalar(prhs[2])+1; // Indexed from 1
			g_alpha = mxGetPr(prhs[3]);
			g_eta = mxGetPr(prhs[4]);
			g_gamma = mxGetPr(prhs[5]);
		}
		else
		{
			printf("Usage : hlda2004(filename,depth,nterms,alpha,eta,gamma,seed(optional))\n");
			return;
		}

		// Seed
		if (nrhs < 7)
			srand((unsigned) time(NULL));
		else
			srand(mxGetScalar(prhs[6]));
		urand();



		mexPrintf("Reading...\n");
		mxGetString(prhs[0],filename,mxGetN(prhs[0])+10);
		docs =  readfile(filename);
		plhs[1] = mxCreateDoubleMatrix(PARAM_MAXITER,1,mxREAL);
		scorearray = mxGetPr(plhs[1]);

		// Multiple initialization suggested in Blei's paper
		bestscore = 0;
		bestdocs = NULL;
		for(iter=0;iter<PARAM_INITITER;iter++)
		{
			topics=reset(docs);
			//DundarSinglePathInitial(docs,topics);
			//halidMultiInitial(docs,topics);
			halidInitial(docs,topics);
			//BleiInitial(docs,topics);
			score = etaScore(topics) + crpScore(topics) + alphaScore(docs);
			mexPrintf("Init iter : %d  Score : %.2f\n",iter,score);
			mexEvalString("drawnow");
			
			if (score > bestscore || iter==0)
			{
				if (iter!=0)
				{
				while(besttopics->count>0)
					removeTopic(besttopics,besttopics->head);
				mxFree(besttopics->head); // Allocated as batch
				}
				besttopics = copyTopics(topics);

				if (bestdocs!=NULL)
					freeDocs(bestdocs);
				bestdocs   = copyDocs(docs);
				bestscore = score;
				plhs[0] = getParentMatrix(topics);
				plhs[2] = getPaths(docs);
			}
			

			
			while(topics->count>0)
				removeTopic(topics,topics->head);
			

			
		}

		topics = besttopics;
		docs = bestdocs;

		// Sample
		lasttime = 0; // print timer
		for(iter=0;iter<PARAM_MAXITER;iter++)
		{
			for (i=0;i<g_ndocs;i++)
			{
			samplePath(&docs[i],topics,0);	
			sampleLevels(&docs[i],0);
			}


			score = etaScore(topics) + crpScore(topics) + alphaScore(docs);
			scorearray[iter] = score;
			if (bestscore < score )
			{
				bestscore = score;
				plhs[0] = getParentMatrix(topics);
				plhs[2] = getPaths(docs);
				plhs[3] = getTopics(topics);
			}
			if ((iter % 5) == 0 && ((time(NULL)-lasttime)>1)) // 1 second between prints
			{
			mexPrintf("Iter:%d nTopics:%d Gibbs Score:%.2f\n",iter,topics->count,score);
			mexEvalString("drawnow");
			lasttime = time(NULL);
			}
			

/*			if (iter%100==0)
			{
				sprintf(filename,"plots/trees/%d/tree%d.txt",topiter,iter);
				file = fopen(filename,"w");

				fprintf(file,"%f\n",score);

				for (i=0,temp=topics->head;i<topics->count;temp=temp->next,i++)
					temp->id=i;
				fprintf(file,"0\t");
				for (i=1,temp=topics->head->next;i<topics->count;temp=temp->next,i++)
					fprintf(file,"%d\t",temp->parent->id);

				fprintf(file,"\n");
				for (i=0,temp=topics->head;i<topics->count;temp=temp->next,i++)
					fprintf(file,"%d\t",temp->ndocs);

				fclose(file);
			}*/

		}

		mexPrintf("Best Gibbs Score:%.2f\n",bestscore);

		

		freeDocs(docs);
		while(topics->count>0)
			removeTopic(topics,topics->head);
		

}
