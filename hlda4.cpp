#include "mex.h"
#include "matrix.h"
#include <string.h>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <list>
#include <algorithm>  
#include "lightspeed\util.h"

using namespace std;

#define PARAM_MAXITER 5000
#define PARAM_INITITER 50

//Global variables
int g_nterms;
int g_depth;
int inmh =0;

class Utils
{
	// Dr. Al Hasan's Algorithm lecture
public:
	template <class T>
	static void inPlacePerm(vector<T>& data)
	{
		int i,r,size;
		T temp;
		size = data.size();
		for (i=0;i<size;i++)
		{
			r = i + ( rand() % (size-i) );
			temp = data[i];
			data[i] = data[r];
			data[r] = temp;
		}
	}

	static vector<int> randomPerm(int size)
	{
		vector<int> perm;
		int i;
		for(i=0;i<size;i++)
		{
			perm.push_back(i);
		}

		inPlacePerm(perm);
		return perm;
	}

	static double urand()
	{
		return (rand()+0.000001) / ((double) (RAND_MAX+1.0));
	}

	template < class T>
	static void printVector(vector<T> vec,string modifier)
	{
		for( vector<T>::const_iterator i = vec.begin(); i != vec.end(); ++i)
			printf(modifier.c_str(),*i);

		printf("\n");
	}

	template < class T>
	static void printVector(vector<T> vec)
	{
		printVector(vec,"\t%d");
	}

	template < class T>
	static void printVector(string name,vector<T> vec)
	{
		printf(name.c_str());
		printVector(vec,"\t%d");
	}



	// ln(Gamma(x+eta)/Gamma(y+eta))
	static double singlelogGammaRatio(int x,int y,double eta)
	{
		int i;
		double result=0;
		if (x-y<6) // Just for efficiency
		{
			for (i=1;i<=(x-y);i++)
			{
			result += log(x-i+eta);
			}
		}
		else result = gammaln(x+eta)-gammaln(y+eta);

		return result;
	}

	template < class T>
	static vector<T> converttoVector(T* data,int size)
	{
		vector<T> res;
		int i;
		for(i=0;i<size;i++)
		{
			res.push_back(data[i]);
		}
		return res;
	}
};


class Doc;
// Default Copy Constructor is sufficient
class Topic
{
public:
	vector<int> words; // Frequency matrix of topic
	int ndocs;
	int nwords; // Total number of words
	int id;		// Topic ID
	double loglikelihood; // Gamma ratio

	Topic()
	{
		ndocs=0;
		nwords=0;
		id=0;
		loglikelihood=0;
		words.resize(g_nterms);
	}

	Topic(const Topic& t): 
	words(t.words),ndocs(t.ndocs),nwords(t.nwords),id(t.id),loglikelihood(t.loglikelihood)
	{
	}

	void print()
	{
		printf("TopicID:\t%d\n",id);
		printf("Nwords:\t%d\n",nwords);
		printf("Ndocs:\t%d\n",ndocs);
	}

};

class hLDATopic : public Topic
{
	public:
	list<hLDATopic>::iterator parent; // Parent node
	list<hLDATopic>::iterator copy; // Facilate copy
	int level;		   // Level in the tree
	int restriction;

	hLDATopic():Topic()
	{
		level = 0;
		restriction =1;
	}

	hLDATopic(const hLDATopic& topic) : Topic(topic)
	{
		level = topic.level;
		parent = topic.parent;
		restriction = topic.restriction;
	}

	void addWords(Doc& d); // Prototype
	void removeWords(Doc& d);


};


class Doc
{
public:
	vector<int> words;
	vector<int> z;			  // Topic assignments path+LDA
	vector<vector<int>> bow;  // Bag of Words, levels , full matrix
	vector<int> indexset;	  // Index of words , sparse index
	int ndiffwords;			  // Number of different words in level
	int length;				  // Document length
	vector<list<hLDATopic>::iterator> path;	  // Document path in the tree
	vector<int> topicsum;     // first g_depth is path topics rest of them is other topics
	int id;

	Doc(int length,int id)
	{
		init();
		words.resize(length);
		z.resize(length);
		this->length = length;
		this->id     = id;
	}

	Doc(string s,int id)
	{
		init();
		readSparse(s);
		this->id = id;
	}

	Doc(string s)
	{
		init();
		readSparse(s);
	}

	// Copy constructor use copied topics , 
	// so topics should be copied first to avoid possible dangling
	// pointer. However in other condition doc copies point same topics
	Doc(const Doc& d) : words(d.words) , z(d.z) , bow[d.bow) , topicsum	(d.topicsum),
		indexset(d.indexset) ,ndiffwords(d.ndiffwords),length(d.length),id(d.id),path(d.path)
	{
		int i;
		//path.resize(g_depth);
	}

	// #items id:frequency format
	void readSparse(string s)
	{
		int wordid,wordcount,i;
		char colon;

		stringstream ss(s);
		ss >> ndiffwords;
		while ((ss >> wordid >> colon >> wordcount))
		{
			length += wordcount;
			//wordid--; // Matlab to C index
			indexset.push_back(wordid); 
			words.reserve(length);
			z.reserve(length);
			for (i=0;i<wordcount;i++)
			{
				words.push_back(wordid);
				z.push_back(0);
			}
		}
		Utils::inPlacePerm(words); // Purmute words
	}

	void init()
	{
		int i;
		ndiffwords=0;
		length = 0;
		id=0;
		path.resize(g_depth);
		topicsum.resize(g_depth);

		bow.resize(g_depth); //2D vectors 
		for (i=0;i<g_depth;i++)
			bow[i].resize(g_nterms);
	}

	// Direct frequency matrix
	void readCSV(string s)
	{
		// TODO: Not implemented yet
	}

	// Read from Matlab Mat array, frequency matrix
	void readMat()
	{
		// TODO: Not implemented yet
	}


	void print()
	{
		int i;

		printf("Docid:\t%d\n",id);
		printf("Length:\t%d\n",length);
		printf("Diffwords:\t%d\n",ndiffwords);
		Utils::printVector("Words:",words);

		printf("Path:");
		//Utils::printVector(path);
		for (i=0;i<g_depth;i++)
		{
			printf("\t%d",path[i]->id);
		}
		printf("\n\n");
		// TODO: Add more info to print
	}

	void printLevels()
	{
		int i,j,sum;
		for(j=0;j<g_depth;j++)
		{
			sum=0;
			for (i=0;i<g_nterms;i++)
				sum += bow[j][i];
			printf("%d\t",sum);
		}
		printf("\n");
		Utils::printVector(topicsum);
	}
};



class hLDATree
{
public:
	list<hLDATopic> topics; // hLDA topics 

	void addToPath(Doc& d,list<hLDATopic>::iterator topic)
	{
		int i,j;
		list<hLDATopic>::iterator temp,prev;
		for (i=topic->level,temp=topic;i>=0;i--)
		{
			d.path[temp->level] = temp;
			temp->ndocs++;
			temp->addWords(d);
			temp = temp->parent;
		}
		// Create new topics  if selected topic is internal node
		for (i=topic->level+1,prev=topic;i<g_depth;i++)
		{
			hLDATopic t;
			t.level = i;
			t.ndocs++  ; // First doc
			t.parent = prev;
			t.addWords(d);
			topics.push_back(t);
			d.path[t.level] = topics.end();
			d.path[t.level]--; //Previous one
			prev = d.path[t.level];
		}
	}

	void removeFromPath(Doc& d)
	{
		int i,j;
		for (i=g_depth-1;i>=0;i--) // From leaf to parent
		{
			d.path[i]->ndocs--;
			d.path[i]->removeWords(d);
			//Remove unnecessary nodes except root
			if (d.path[i]->ndocs == 0 && i!=0) 
			{
				topics.erase(d.path[i]); // It is no longer valid iterator
			}
		}
	
	}

	void print()
	{
		list<hLDATopic>::iterator it;
		int i;
		printf("\nTopics:");
		for (i=0,it=topics.begin();it!=topics.end();it++,i++)
			it->id = i;
		for (it=topics.begin();it!=topics.end();it++)
			printf("\t%d",it->id);
		printf("\n\n");
		for (it=topics.begin();it!=topics.end();it++)
			printf("\t%d",it->parent->id);
		printf("\n\n");
		for (it=topics.begin();it!=topics.end();it++)
			printf("\t%d",it->ndocs);
	}

	void relax()
	{
		list<hLDATopic>::iterator temp;
		for(temp=topics.begin();temp!=topics.end();temp++)
		{
			temp->restriction = 1;
		}
		
	}
};





class State
{
	public:
	int ndocs;				 // Number of documents in the corpus 
	vector<double> eta;		 // Topic Dirichlet Param
	vector<double> gamma;	 // nCRP param
	vector<double> alphal;	 // Local LDA document Dirichlet param
	hLDATree	   tree;	 // hLDA Tree
	vector<Doc>    docs;     // Documents

	State(vector<double> alphal,vector<double> eta,vector<double> gamma)
	{
		ndocs=0;
		this->alphal = alphal;
		this->eta    = eta;
		this->gamma  = gamma;
		hLDATopic root;
		tree.topics.push_back(root);
		tree.topics.begin()->parent = tree.topics.begin();
	}

	State(State &s)
	{
		copy(s);
	}

	State State::operator=(State& s)
	{
		copy(s);
		return *this;
	}

	void copy(State &s)
	{
	    eta= s.eta;
		ndocs=s.ndocs;
		gamma=s.gamma;
		tree=s.tree;
		docs=s.docs;
		alphal=s.alphal;

		int i,j;
		list<hLDATopic>::iterator iter,itmez;
		
		for(i=0,
			iter=s.tree.topics.begin(),
			itmez=tree.topics.begin()
			;i<s.tree.topics.size();
			i++,iter++,itmez++)//Copy iterators
			{
			iter->copy = itmez;
			itmez->parent = iter->parent->copy;
			}

		for(i=0;i<ndocs;i++)//Copy iterators
		{
			for (j=0;j<g_depth;j++)
			{
				docs[i].path[j]	= docs[i].path[j]->copy;
			}
		}
	}
	
	
	// ntopics x (g_nterms+2)
	mxArray* getTreeTopics()
	{
		mxArray* parents=mxCreateDoubleMatrix(tree.topics.size(),g_nterms+2,mxREAL);
		list<hLDATopic>::iterator it;
		int i,j;
		double*  p = mxGetPr(parents);
		for(i=0,it=tree.topics.begin();it!=tree.topics.end();it++,i++)
		{
			it->id = i+1; // new ids
		}
		p[0] = 0;
		for(i=0,it=tree.topics.begin();it!=tree.topics.end();it++,i++)
		{
			if (i!=0)
				p[i] = it->parent->id;
			p[i+tree.topics.size()] = it->ndocs;
			for (j=0;j<g_nterms;j++)
				p[i+(tree.topics.size())*(j+2)] = it->words[j];
		}
		return parents;
	}

};



//**  Member Functions

void hLDATopic::addWords(Doc& d)
{
		int i,j,level;
		level = this->level;
		for(j=0;j<d.indexset.size();j++)
		{
			words[d.indexset[j]] += d.bow[level][d.indexset[j]];
		}
		nwords += d.topicsum[level];
}


void hLDATopic::removeWords(Doc& d)
{
		int i,j,level;
		level = this->level;
		for(j=0;j<d.indexset.size();j++)
		{
			words[d.indexset[j]] -= d.bow[level][d.indexset[j]];
		}
		nwords -= d.topicsum[level];
}


class Sampler
{
public:
	State state;


	Sampler(vector<double> alphal,vector<double> eta,vector<double> gamma) : state(State(alphal,eta,gamma))
	{
	}


	double sampleLevels(Doc& d,int init=0)
	{
		vector<double> topicpresence(g_depth),topiclikelihood(g_depth),prob(g_depth);
		double denominator,arand,loglikelihood;
		int i,j,word,topic;
		
		loglikelihood=0;


		for(i=0;i<d.length;i++)
		{
			word = d.words[i];
			topic = d.z[i];

			// Remove word from topic
			if (init!=1)
			{
				d.path[topic]->words[word]--;
				d.path[topic]->nwords--;
				d.bow[topic][word]--;
				d.topicsum[topic]--;

			}

			//HLDA topics
			denominator = 0;
			for (j=0;j<g_depth;j++)
			{
				topicpresence[j] = d.topicsum[j] + state.alphal[j];
				topiclikelihood[j] = ( d.path[j]->words[word] + state.eta[j]) / ( d.path[j]->nwords + g_nterms*state.eta[j]) ;
				prob[j] = topicpresence[j] * topiclikelihood[j];
				denominator+=prob[j];
			}


			//Sample
			arand = Utils::urand() * denominator;
			for(j=0;j<g_depth;j++)
			{
				if 	(arand <= prob[j])
				{
					topic = j;
					break;
				}
				else
					arand = arand - prob[j];
			}

			loglikelihood += log(prob[topic]/denominator);
			d.topicsum[topic]++;
			d.z[i] = topic;
			d.path[topic]->words[word]++;
			d.path[topic]->nwords++;
			d.bow[topic][word]++;
		}
		return loglikelihood;
	}


	double imitateLevels(Doc& d,Doc& reference)
	{
		vector<double> topicpresence(g_depth),topiclikelihood(g_depth),prob(g_depth);
		double denominator,arand,loglikelihood;
		int i,j,s,word,topic;
		
		loglikelihood=0;

		vector<int> perm = Utils::randomPerm(d.length);

		for(s=0;s<d.length;s++)
		{
			i = perm[s];
			word = d.words[i];
			topic = d.z[i];

			d.path[topic]->words[word]--;
			d.path[topic]->nwords--;
			d.bow[topic][word]--;
			d.topicsum[topic]--;

			//HLDA topics
			denominator = 0;
			for (j=0;j<g_depth;j++)
			{
				topicpresence[j] = d.topicsum[j] + state.alphal[j];
				topiclikelihood[j] = ( d.path[j]->words[word] + state.eta[j]) / ( d.path[j]->nwords + g_nterms*state.eta[j]) ;
				prob[j] = topicpresence[j] * topiclikelihood[j];
				denominator+=prob[j];
			}


			//Imitate
			topic = reference.z[i];

			loglikelihood += log(prob[topic]/denominator);
			d.topicsum[topic]++;
			d.z[i] = topic;
			d.path[topic]->words[word]++;
			d.path[topic]->nwords++;
			d.bow[topic][word]++;
		}
		return loglikelihood;
	}

	//Empty topic gamma ratio , cumulative
	vector<double> gammaRatio(Doc& d)
	{
		int w,i,l;
		vector <double> cumlogsum(g_depth);

		cumlogsum[g_depth-1] = 0; 
		for(l=g_depth-1;l>0;l--) // From reverse , root has to add all empty topics
		{
			cumlogsum[l-1] = cumlogsum[l];
			for (i=0;i<d.ndiffwords;i++)
			{
				w = d.indexset[i];
				if (d.bow[l][w]==0) continue;
				cumlogsum[l-1] += Utils::singlelogGammaRatio(d.bow[l][w],0,state.eta[l]);
			}
			cumlogsum[l-1] -= Utils::singlelogGammaRatio(d.topicsum[l],0,state.eta[l]*g_nterms);
		}
		
		return cumlogsum;
	}

	double gammaRatio(Doc& d,list<hLDATopic>::iterator t)
	{
		int w,i,l=t->level;
		double logsum=0;
		for (i=0;i<d.ndiffwords;i++)
		{
			w = d.indexset[i];
			if (d.bow[l][w]==0) continue;
			logsum += Utils::singlelogGammaRatio(d.bow[l][w]+t->words[w],t->words[w],state.eta[l]);
		}
		logsum -= Utils::singlelogGammaRatio(d.topicsum[l]+t->nwords,t->nwords,state.eta[l]*g_nterms);
		
		return logsum;
	}


	double samplePath(Doc& d,int init,list<hLDATopic>::iterator merge)
	{
		int i;
		int size = state.tree.topics.size();
		list<hLDATopic>::iterator temp;
		for(temp=state.tree.topics.begin();temp!=state.tree.topics.end();temp++)
		{
			if (temp==merge || temp->parent == merge) // 3 level
			{ 
				temp->restriction = 1;
			}
			else
			{
				temp->restriction = 0; 
			}
		}
		return samplePath(d,init);
	}



	double samplePath(Doc& d,int init,list<hLDATopic>::iterator split1,list<hLDATopic>::iterator split2)
	{
		int i;
		int size = state.tree.topics.size();
		list<hLDATopic>::iterator temp;
		for(temp=state.tree.topics.begin();temp!=state.tree.topics.end();temp++)
		{
			if (temp==split1 || temp->parent == split1 || temp==split2 || temp->parent == split2) // 3 level
			{ 
				temp->restriction = 1;
			}
			else
			{
				temp->restriction = 0; 
			}
		}
		return samplePath(d,init);
	}

	


	double samplePath(Doc& d,int init)
	{

		vector<double> empty,ncrp,likelihood,gammaratio,probs; //All in logs
		double maxprob,denominator,arand,pval;
		vector<list<hLDATopic>::iterator>::iterator topic;
		list<hLDATopic>::iterator candidate,temp;
		int i,j;

	

		if (!init)
		{
			state.tree.removeFromPath(d);
		}

		likelihood.resize(state.tree.topics.size());
		ncrp.resize(state.tree.topics.size());
		probs.resize(state.tree.topics.size());

		empty = gammaRatio(d);

		for(candidate=state.tree.topics.begin();candidate!=state.tree.topics.end();candidate++)
		{
			candidate->loglikelihood = gammaRatio(d,candidate);
		}

		// Collect probabilities
		for(i=0,candidate=state.tree.topics.begin();candidate!=state.tree.topics.end();candidate++,i++)
		{
			ncrp[i] = 0;
			likelihood[i] = 0;

			for(j=candidate->level,temp=candidate;j>0;temp=temp->parent,j--)
			{
				// Use its parents ncrp distribution
				ncrp[i] += log(temp->ndocs / (temp->parent->ndocs + state.gamma[temp->level-1] ));
			}

			// Open new path if internal node selected
			if(candidate->level < g_depth-1)
				ncrp[i] += log( state.gamma[candidate->level] / (candidate->ndocs + state.gamma[candidate->level] ));

			for(j=candidate->level,temp=candidate;j>=0;temp=temp->parent,j--)
			{
				likelihood[i] += temp->loglikelihood;
			}

			likelihood[i] += empty[candidate->level];
			probs[i] = ncrp[i] + likelihood[i];
		}

		denominator = 0;
		maxprob = probs.back(); // For split merge
		arand   = Utils::urand();

		/*if (inmh==1)
		{
		printf("%f\n",arand);
		}*/

		for(temp=state.tree.topics.begin(),i=0;temp!=state.tree.topics.end();temp++,i++)
		{
			if (maxprob<probs[i]&&temp->restriction==1)
				maxprob = probs[i]; //Restrict some specified nodes
		}

		for(temp=state.tree.topics.begin(),i=0;temp!=state.tree.topics.end();temp++,i++)
		{
			probs[i] -= maxprob;
			probs[i] = exp(probs[i])*temp->restriction; // Convert to actual prob
			denominator += probs[i];
		}

	
		arand = arand*denominator;

/*		if (inmh==1)
		{
		printf("\n");
		for(i=0,candidate=state.tree.topics.begin();candidate!=state.tree.topics.end();candidate++,i++)
		{
			printf("%d ",candidate->restriction);
		}
		printf("\n");
		}
		*/
		for(i=0,candidate=state.tree.topics.begin();candidate!=state.tree.topics.end();candidate++,i++)
		{
			if (probs[i]>=arand)
			{
				temp = candidate;
				break;
			}
			else
				arand -= probs[i];
		}


		int level = temp->level;

		state.tree.addToPath(d,temp);
		

		return log(probs[i]/denominator);
	}

	double imitatePath(int d,State& s,State& reference,vector<list<hLDATopic>::iterator>& map1,vector<list<hLDATopic>::iterator>& map2)
	{

		vector<double> empty,ncrp,likelihood,gammaratio,probs; //All in logs
		double maxprob,denominator,arand,pval;
		vector<list<hLDATopic>::iterator>::iterator topic;
		list<hLDATopic>::iterator candidate,temp,selected;
		int i,j,inlist;

		s.tree.removeFromPath(s.docs[d]);


		likelihood.resize(s.tree.topics.size());
		ncrp.resize(s.tree.topics.size());
		probs.resize(s.tree.topics.size());

		empty = gammaRatio(s.docs[d]);

		for(candidate=s.tree.topics.begin();candidate!=s.tree.topics.end();candidate++)
		{
			candidate->loglikelihood = gammaRatio(s.docs[d],candidate);
			candidate->restriction = 0;
		}

		// Collect probabilities
		for(i=0,candidate=s.tree.topics.begin();candidate!=s.tree.topics.end();candidate++,i++)
		{
			ncrp[i] = 0;
			likelihood[i] = 0;

			for(j=candidate->level,temp=candidate;j>0;temp=temp->parent,j--)
			{
				// Use its parents ncrp distribution
				ncrp[i] += log(temp->ndocs / (temp->parent->ndocs + s.gamma[temp->level-1] ));
			}

			// Open new path if internal node selected
			if(candidate->level < g_depth-1)
				ncrp[i] += log( s.gamma[candidate->level] / (candidate->ndocs + s.gamma[candidate->level] ));

			for(j=candidate->level,temp=candidate;j>=0;temp=temp->parent,j--)
			{
				likelihood[i] += temp->loglikelihood;
			}

			likelihood[i] += empty[candidate->level];
			probs[i] = ncrp[i] + likelihood[i];
		}

		denominator = 0;
		maxprob = probs.back();
		arand   = Utils::urand();



		for(i=0,inlist=0;i<map1.size();i++)
			map1[i]->restriction = 1;

		for(i=0,temp=s.tree.topics.begin();i<s.tree.topics.size();temp++,i++)
		{
			if (maxprob<probs[i]&&temp->restriction)
				maxprob = probs[i]; //Restrict some specified nodes
		}

		for(temp=s.tree.topics.begin(),i=0;temp!=s.tree.topics.end();temp++,i++)
		{
			probs[i] -= maxprob;
			probs[i] = exp(probs[i])*temp->restriction; // Convert to actual prob
			denominator += probs[i];
		}

		// Imitate Sampling
		// I need index of topic in current state 
		// based on topic index in refence 
		// We can obtain topic index in refence by for loop
		// But how to match topic index in current state. 
		for(i=0,inlist=0;i<map1.size();i++)
		{
			if (reference.docs[d].path[2]==map2[i])
			{
				inlist = 2;
				break;
			}
		}
		if (inlist == 2)
		{
			s.tree.addToPath(s.docs[d],map1[i]);	
			selected=map1[i];
		}
		else
		{
			for(i=0,inlist=0;i<map1.size();i++)
			{
				if (reference.docs[d].path[1]==map2[i])
				{
					inlist = 1;
					break;
				}
			
			}
			selected=map1[i];
			s.tree.addToPath(s.docs[d],map1[i]);
			map1.push_back(s.docs[d].path[2]);
			map2.push_back(reference.docs[d].path[2]);
		}
		


		s.tree.relax();
		for(temp=s.tree.topics.begin(),i=0;temp!=s.tree.topics.end();temp++,i++)
		{
			if (temp==selected)
				return log(probs[i]/denominator);
		}
		return 1;
	}

	double etaScore()
	{
		int i;
		double logsum;
		list<hLDATopic>::iterator topic;

		logsum=0;
		for(topic=state.tree.topics.begin();topic!=state.tree.topics.end();topic++)
		{
			for (i=0;i<g_nterms;i++)
			{
				if(topic->words[i]==0)	continue; // There could be better optimization here
			logsum += Utils::singlelogGammaRatio(topic->words[i],0,state.eta[topic->level]);
			}
			logsum -= Utils::singlelogGammaRatio(topic->nwords,0,g_nterms*state.eta[topic->level]);
		}

		return logsum;
	}

	double alphaScore()
	{
		int i,j;
		double logsum;
		double alphasum,sumgammalnalpha;

		alphasum = 0;
		sumgammalnalpha = 0;
		logsum = 0;

		for(j=0;j<g_depth;j++)
		{
			alphasum += state.alphal[j];
		}

		for(i=0;i<state.docs.size();i++)
		{
			logsum -= Utils::singlelogGammaRatio(state.docs[i].length,0,alphasum);
			for(j=0;j<g_depth;j++){
				logsum += Utils::singlelogGammaRatio(state.docs[i].topicsum[j],0,state.alphal[j]);
			}
		}
		return logsum;
	}

	double crpScore()
	{
		double logsum= - Utils::singlelogGammaRatio(state.tree.topics.begin()->ndocs,0,state.gamma[0]); //Root
		list<hLDATopic>::iterator topic;
		for(topic=state.tree.topics.begin(),topic++;topic!=state.tree.topics.end();topic++) //Start from 2nd
		{
			if(topic->level<g_depth-1) // Non leaf
			{
				logsum -= Utils::singlelogGammaRatio(topic->ndocs,0,state.gamma[topic->level]);
			}
			logsum += log(state.gamma[topic->parent->level]) + gammaln(topic->ndocs);
		}
		return logsum;
	}



	double score()
	{
		return etaScore() + crpScore() + alphaScore();
	}

	void halidInitial()
	{
		int i,j;
		for (i=0;i<state.ndocs;i++)
		{
			for(j=0;j<state.docs[i].length;j++)
			{
				state.docs[i].z[j] = (rand() % (g_depth)); // Level starts from 0
				state.docs[i].bow[ state.docs[i].z[j]][state.docs[i].words[j]]++;
				state.docs[i].topicsum[state.docs[i].z[j]]++;
			}
		}

		for(i=0;i<state.docs.size();i++)
		{
			
			samplePath(state.docs[i],1);
			sampleLevels(state.docs[i],0);
		}

	}

	void sample()
	{
		int i,ndocs;
		ndocs = state.docs.size();


			for (i=0;i<ndocs;i++)
			{
				sampleLevels(state.docs[i],0);
				samplePath(state.docs[i],0);
			}
			printf("Gibbs->:%d\n",state.tree.topics.size());
			mh();
			printf("MH->:%d\n",state.tree.topics.size());
	}

	void mh()
	{
		int i,j,k,d;
		double loglikelihood,yedek;

		i = rand() % state.ndocs;
		j = rand() % state.ndocs;

		
		State copy=state,launch(state);

		list<hLDATopic>::iterator root,clust1,clust2,merge;

		vector<int> S;
		vector<int> restrict;

		root   = state.tree.topics.begin();
		clust1 = state.docs[i].path[1];
		clust2 = state.docs[j].path[1];


		double oldscore = score();

		if (i!=j)
		{
			inmh = 1;
			if (clust1==clust2)
			{
				printf("Split\n");
				for(k=0;k<state.ndocs;k++)
				{
					if (k!=j&&k!=i &&  ( state.docs[k].path[1] == clust1 ))
					{
						S.push_back(k);
					}
				}
				Utils::inPlacePerm(S);

				//printf("NTOPICS:%d\n",state.tree.topics.size());
				
				state.tree.removeFromPath(state.docs[i]);
				state.tree.addToPath(state.docs[i],root);
				loglikelihood = sampleLevels(state.docs[i],0);

				state.tree.removeFromPath(state.docs[j]);
				//printf("NTOPICS:%d\n",state.tree.topics.size());
				state.tree.addToPath(state.docs[j],root);
				loglikelihood += sampleLevels(state.docs[j],0);

				/*for(k=0;k<S.size();k++)
				{
					d=S[k];
					samplePath(state.docs[d],0,state.docs[j].path[1],state.docs[i].path[1]);
					sampleLevels(state.docs[d],0);
				}

				launch = state;*/


				for(k=0;k<S.size();k++)
				{
					d=S[k];
					loglikelihood += samplePath(state.docs[d],0,state.docs[j].path[1],state.docs[i].path[1]);
					//printf("Path:%f\n",loglikelihood);
					loglikelihood += sampleLevels(state.docs[d],0);
					//printf("Level:%f\n",loglikelihood);
					/* state.tree.removeFromPath(state.docs[d]);
						state.tree.addToPath(state.docs[d],merge);*/
				}
				//launch.tree.print();
				for (int at=0;at<2;at++)
				{
				sampleLevels(state.docs[i],0);
				sampleLevels(state.docs[j],0);
				for(k=0;k<S.size();k++)
				{
					d=S[k];
					sampleLevels(state.docs[d],0);
				}
				}

				//printf("NTOPICS:%d\n",state.tree.topics.size());

				loglikelihood = -loglikelihood; // denominator
				//Ad indices of the i,j 
				S.push_back(i);
				S.push_back(j);
				loglikelihood += reverseSplitLikelihood(copy,state,S);
				


				
				//copy.tree.print();
				//state.tree.print();

				if (exp(score()-oldscore) < Utils::urand())
				{
					state = copy;
				}
				else
				{

					printf("\nAccepted\n");
					state.tree.relax();
				}

			}
			else
			{
				printf("Merge\n");
				for(k=0;k<state.ndocs;k++)
				{
					if (k!=j&&k!=i &&  ( state.docs[k].path[1] == clust1 || state.docs[k].path[1] == clust2))
					{
						S.push_back(k);
					}
				}
				Utils::inPlacePerm(S);
				
				//printf("NTOPICS:%d\n",state.tree.topics.size());

				state.tree.removeFromPath(state.docs[i]);
				state.tree.addToPath(state.docs[i],root);
				merge = state.docs[i].path[1];
				//sampleLevels(state.docs[i],0);
				loglikelihood = sampleLevels(state.docs[i],0);

				
				//printf("NTOPICS:%d\n",state.tree.topics.size());
				//samplePath(state.docs[j],0,merge);
				//sampleLevels(state.docs[j],0);
				loglikelihood += samplePath(state.docs[j],0,merge);
				
				loglikelihood += sampleLevels(state.docs[j],0);


				/*for(k=0;k<S.size();k++)
				{
					d=S[k];
					samplePath(state.docs[d],0,merge);
					sampleLevels(state.docs[d],0);
				}

				launch = state;
				*/

				

				for(k=0;k<S.size();k++)
				{
					d=S[k];
					loglikelihood += samplePath(state.docs[d],0,merge);
					//printf("Path:%f\n",loglikelihood);

					loglikelihood += sampleLevels(state.docs[d],0);
					//printf("Level:%f\n",loglikelihood);
					/* state.tree.removeFromPath(state.docs[d]);
						state.tree.addToPath(state.docs[d],merge);*/
				}
				

				for (int at=0;at<20;at++)
				{
				sampleLevels(state.docs[i],0);
				sampleLevels(state.docs[j],0);
				for(k=0;k<S.size();k++)
				{
					d=S[k];
					sampleLevels(state.docs[d],0);
				}
				}
				//launch.tree.print();
				//state.tree.print();
				//printf("NTOPICS:%d\n",state.tree.topics.size());

				loglikelihood = -loglikelihood; // denominator
				yedek = loglikelihood;
				//Ad indices of the i,j 
				S.push_back(i);
				S.push_back(j);
				loglikelihood += reverseMergeLikelihood(copy,state,S);


				
				
				state.tree.relax();
				loglikelihood -= oldscore;
				loglikelihood += score();
				//printf("Accept:%f\n",(loglikelihood));

				//copy.tree.print();
				//state.tree.print();


				if (exp(score()-oldscore) < Utils::urand())
				{
					state = copy;
				}
				else
				{
					copy.tree.print();
					state.tree.print();
					printf("\nAccepted\n");
					state.tree.relax();
				}
			}
		}
	inmh =0;
	}

	double reverseMergeLikelihood(State old,State s,vector<int> docs)
	{
		double loglikelihood = 0;
		int i,j,k,d,level,doc1,doc2;
		vector<double> prior(g_depth),likelihood(g_depth);
		vector<list<hLDATopic>::iterator> map1,map2;

		
		loglikelihood = 0;

		// Create split configuration
		doc2 = docs.back();
		docs.pop_back();
		doc1 = docs.back();
		docs.pop_back();

		s.tree.removeFromPath(s.docs[doc1]);
		s.tree.addToPath(s.docs[doc1],s.tree.topics.begin());
		map1.push_back(s.docs[doc1].path[1]);
		map1.push_back(s.docs[doc1].path[2]);
		map2.push_back(old.docs[doc1].path[1]);
		map2.push_back(old.docs[doc1].path[2]);
		loglikelihood = imitateLevels(s.docs[doc1],old.docs[doc1]);

		s.tree.removeFromPath(s.docs[doc2]);
		s.tree.addToPath(s.docs[doc2],s.tree.topics.begin());
		map1.push_back(s.docs[doc2].path[1]);
		map1.push_back(s.docs[doc2].path[2]);
		map2.push_back(old.docs[doc2].path[1]);
		map2.push_back(old.docs[doc2].path[2]);
		loglikelihood += imitateLevels(s.docs[doc2],old.docs[doc2]);

		Utils::inPlacePerm(docs);

		for (i=0;i<docs.size();i++)
		{
			d = docs[i];
			loglikelihood += imitatePath(d,s,old,map1,map2);
			loglikelihood += imitateLevels(s.docs[d],old.docs[d]);
			
			// Imitate path sampling

		}
		return loglikelihood;
	}

	double reverseSplitLikelihood(State old,State s,vector<int> docs)
	{
		double loglikelihood = 0;
		int i,j,k,d,level,doc1,doc2;
		vector<double> prior(g_depth),likelihood(g_depth);
		vector<list<hLDATopic>::iterator> map1,map2;

		
		loglikelihood = 0;

		// Create split configuration
		doc2 = docs.back();
		docs.pop_back();
		doc1 = docs.back();
		docs.pop_back();

		s.tree.removeFromPath(s.docs[doc1]);
		s.tree.addToPath(s.docs[doc1],s.tree.topics.begin());
		map1.push_back(s.docs[doc1].path[1]);
		map1.push_back(s.docs[doc1].path[2]);
		map2.push_back(old.docs[doc1].path[1]);
		map2.push_back(old.docs[doc1].path[2]);
		loglikelihood = imitateLevels(s.docs[doc1],old.docs[doc1]);

		loglikelihood += imitatePath(doc2,s,old,map1,map2);
		loglikelihood += imitateLevels(s.docs[doc2],old.docs[doc2]);

		Utils::inPlacePerm(docs);

		for (i=0;i<docs.size();i++)
		{
			d = docs[i];
			loglikelihood += imitatePath(d,s,old,map1,map2);
			loglikelihood += imitateLevels(s.docs[d],old.docs[d]);
			
			// Imitate path sampling

		}
		return loglikelihood;
	}

	
	void readDocs(string filename)
	{
		ifstream input(filename);
		string line;
		int i;

		while (getline (input,line)) // Read docs
		{
			state.docs.push_back(Doc(line));
			state.ndocs++;
		}

		input.close();
	}
};



void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
                 const mxArray *prhs[])
{
	char* filename;
	double lasttime;
	vector<double> alphal,gamma,eta;
	int iter;

	if (nrhs>5) {
			filename = (char*) mxCalloc((mxGetN(prhs[0])+1),sizeof(char));
			g_depth = mxGetScalar(prhs[1]);
			g_nterms = mxGetScalar(prhs[2]); // Indexed from 0
			alphal = Utils::converttoVector(mxGetPr(prhs[3]),g_depth);
			eta	   = Utils::converttoVector(mxGetPr(prhs[4]),g_depth);
			gamma  = Utils::converttoVector(mxGetPr(prhs[5]),g_depth-1);
	}
	else 
	{
			printf("Usage : hlda2(filename,depth,nterms,alphal,eta,gamma,seed(optional))\n");
			return;
	}

	if (nrhs < 7)
			srand((unsigned) time(NULL));
	else
			srand(mxGetScalar(prhs[7]));
	Utils::urand();





	Sampler s(alphal,eta,gamma);
	printf("Reading...\n");
	mexEvalString("drawnow");
	mxGetString(prhs[0],filename,mxGetN(prhs[0])+1);
	s.readDocs(filename);

	printf("Initialization...\n");
	mexEvalString("drawnow");

	s.halidInitial();

	plhs[0]=s.state.getTreeTopics();
	plhs[1] = mxCreateDoubleMatrix(PARAM_MAXITER,1,mxREAL);
	double gibbs,maxscore = s.score();

	lasttime = time(NULL);
	printf("Sampling...\n");
	mexEvalString("drawnow");
	for(iter=0;iter<PARAM_MAXITER;iter++)
	{
		if (iter%10==0)
		{
			//printf("%d-%d-%.2f\n",iter,s.state.tree.topics.size(),s.score());
			//mexEvalString("drawnow");
		}
		s.sample();
		gibbs = s.score();
		mxGetPr(plhs[1])[iter] = gibbs;

		if ((iter % 5) == 0 && ((time(NULL)-lasttime)>1)) // 1 second between prints
		{
			mexPrintf("Iter:%d nTopics:%d Gibbs Score:%.2f\n",iter,s.state.tree.topics.size(),gibbs);
			mexEvalString("drawnow");
			lasttime = time(NULL);
		}
		
		if (gibbs>maxscore)
		{
			plhs[0]=s.state.getTreeTopics();
			
			maxscore = gibbs;
			printf("Best : %f\n",gibbs);
			mexEvalString("drawnow");
		}
	}


}
