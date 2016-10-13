#include "Doc.h"
#include "Topic.h"

void inPlacePerm(vector<int>& data)
{
	int size = data.size();
	int i, temp, r;
	for (i = 0; i<size; i++)
	{
		r = i + (rand() % (size - i));
		temp = data[i];
		data[i] = data[r];
		data[r] = temp;
	}
}



using namespace std;
Doc::Doc(string line,int id)
{
	
	int wordid, wordcount, i;
	char colon;
	validpath = 0;
	stringstream ss(line);
	ss >> ndiffwords;

	bow = Matrix(depth, nterms);
	bow.zero();
	levelsum.resize(depth);
	emptyscore.resize(depth); //  = Vector(depth);
	path.resize(depth);
	prob.resize(depth);

	length = 0;
	for (int i = 0; i < depth; i++)
	{
		levelsum[i] = 0;
	}

	this->id = id;

	while ((ss >> wordid >> colon >> wordcount))
	{
		length += wordcount;

		//wordid--; // -1 if word index starts from 1

		indexset.push_back(wordid);
		words.resize(length);
		z.reserve(length);
		for (i = 0; i<wordcount; i++)
		{
			words[length - i - 1] = wordid; 
			z.push_back(0);
		}
	}
	inPlacePerm(words);
	/*nheld = floor(length / 5.0);
	heldout = vector<int>(words.end() - nheld, words.end());
	length = length - nheld;
	words.resize(length);
	
	heldz.resize(nheld);
	z.resize(length);*/
	//random_shuffle(words.begin(), words.end());
	
}

Doc::~Doc()
{

}




void Doc::samplezfromPrior()
{
	double alphasum = alpha.sum();
	levelsum.zero();
	bow.zero();
	for (int i = 0; i < length;i++)
	{
		z[i] =  ((alpha) / (alphasum)).sample(); // i is the number of sampled words previously  // rand() % depth; //
		levelsum[z[i]]++;
		bow.data[z[i]*nterms + words[i]]++;
	}
}

void Doc::sampleunif() //Uniform sampling
{
	double alphasum = alpha.sum();
	levelsum.zero();
	bow.zero();
	for (int i = 0; i < length; i++)
	{
		z[i] = rand()%depth; // i is the number of sampled words previously  // rand() % depth; //
		levelsum[z[i]]++;
		bow.data[z[i] * nterms + words[i]]++;
	}
}

/*void Doc::sampleunifheld() //Uniform sampling
{
	double alphasum = alpha.sum();
	for (int i = 0; i < nheld; i++)
	{
		heldz[i] = rand() % depth; // i is the number of sampled words previously  // rand() % depth; //
		levelsum[heldz[i]]++;
		bow.data[heldz[i] * nterms + heldout[i]]++;
		path[heldz[i]]->words[heldout[i]]++;
		path[heldz[i]]->nwords++;
	}
}*/

void Doc::allLevel2()
{
	levelsum.zero();
	bow.zero();
	for (int i = 0; i < length; i++)
	{
		z[i] = 2;
		levelsum[z[i]]++;
		bow.data[z[i] * nterms + words[i]]++;
	}
}



void Doc::emptyScore()
{
	int j, w, l;
	double logsum = 0;
	emptyscore[0] = 0; // Root does not have any impact
	double* bow;
	double etal,etasl,lgetal;
	for (l = 1; l < depth; l++)
	{

		logsum = 0;
		bow		= this->bow.data + l*nterms;
		etal = eta[l]; 
		lgetal = gammatable.get(etal, 0);
		etasl	= etasum[l];
		for (j = 0; j < ndiffwords; j++){
			w = indexset[j];
			if (bow[w] == 0) continue;
			logsum += gammatable.get(etal,bow[w]) - lgetal;
		}
		logsum += gammatable.get(etasl, 0) - gammatable.get(etasl, levelsum.data[l]);
		emptyscore[l] = logsum;
	}
}

double Doc::sampleLevels(int calcprob)
{
	int i, j,l, w;
	double sum=0,gibbs=0,r;

	for (i = 0; i < length; i++)
	{
		w = words[i];
		l = z[i];

		
		path[l]->words[w]--;
		path[l]->nwords--;
		levelsum[l]--;
		
		sum = 0;
		//parallel_for(0, depth, [&](int j)
		for (j = 0; j < depth; j++) // It will be faster if memory access is consecutive for topics as in Topic Model toolbox
		{
		// Most computational part
			prob[j] = ((levelsum[j] + alpha[j]) * (path[j]->words[w] + eta[j])) / (path[j]->nwords + etasum[j]);
			sum += prob[j];
		}//);
		// Heavy computational part
		//l = prob.sample();

		r = urand()*sum; // Mersenne twister in cokus.h is faster
		for (j = 0; j < depth; j++)
		{
			if (r >= prob[j])
				r = r - prob[j];
			else
				break;
		}

		//Performance critical , update stats
		if (l != j)
		{
			bow[l][w]--;
			bow[j][w]++; 
		}
		l = j;

		if (calcprob) // Avoid extra computation , no need to calculate likelihood for gibbs
		gibbs += log(prob[l]/sum);
		path[l]->words[w]++;
		path[l]->nwords++;
		levelsum[l]++;
		
		z[i] = l;
	}
	return gibbs;
}
/*
double Doc::sampleLevelsheld(int calcprob)
{
	int i, j, l, w;
	double sum = 0, gibbs = 0, r;

	for (i = 0; i < nheld; i++)
	{
		w = heldout[i];
		l = heldz[i];


		path[l]->words[w]--;
		path[l]->nwords--;
		levelsum[l]--;

		sum = 0;
		//parallel_for(0, depth, [&](int j)
		for (j = 0; j < depth; j++) // It will be faster if memory access is consecutive for topics as in Topic Model toolbox
		{
			// Most computational part
			prob[j] = ((levelsum[j] + alpha[j]) * (path[j]->words[w] + eta[j])) / (path[j]->nwords + etasum[j]);
			sum += prob[j];
		}//);
		// Heavy computational part
		//l = prob.sample();

		r = urand()*sum; // Mersenne twister in cokus.h is faster
		for (j = 0; j < depth; j++)
		{
			if (r >= prob[j])
				r = r - prob[j];
			else
				break;
		}

		//Performance critical , update stats
		if (l != j)
		{
			bow[l][w]--;
			bow[j][w]++;
		}
		l = j;

		if (calcprob) // Avoid extra computation , no need to calculate likelihood for gibbs
			gibbs += log(prob[l] / sum);
		path[l]->words[w]++;
		path[l]->nwords++;
		levelsum[l]++;

		heldz[i] = l;
	}
	return gibbs;
}*/

double Doc::immitateLevels(vector<Doc>::iterator& d)
{
	Vector prob(depth); // Actually there is not need to reserve this each time
	int i, j, l, w;
	double sum = 0, gibbs = 0;
	for (i = 0; i < length; i++)
	{
		w = words[i];
		l = z[i];

		path[l]->words[w]--;
		path[l]->nwords--;
		levelsum[l]--;
		bow.data[l*nterms + w]--; //Performance critical
		sum = 0;
		for (j = 0; j < depth; j++)
		{
			// Most computational part
			prob[j] = (levelsum[j] + alpha[j]) * ((path[j]->words[w] + eta[j]) / (path[j]->nwords + etasum[j]));
			sum += prob[j];
		}
		// Heavy computational part
		l = d->z[i]; // Just immitation
		gibbs += log(prob[l] / sum);
		path[l]->words[w]++;
		path[l]->nwords++;
		levelsum[l]++;
		bow.data[l*nterms + w]++; //Performance critical
		z[i] = l;
	}
	return gibbs;
}