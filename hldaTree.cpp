#include "hldaTree.h"


hldaTree::hldaTree()
{
	topics.emplace_back(0);
	ntopics = 1;
	nextid = 1;
}


hldaTree::~hldaTree()
{
}


void hldaTree::likelihoods(vector<Doc>::iterator doc)
{
	int j, w, l,c,t,tw,etal;
	double loglikelihood;

	for (auto &t : topics){
		l = t.level;
		loglikelihood = 0;

		if (l == 0) { // Root common in all paths
			t.loglikelihood = 0;  continue;
		}

		etal = eta[l];
		double* topic = t.words.data;
		double* bow = doc->bow.data + l*nterms;
		double* table = gammatable[eta[l]].data;
		for (j = 0; j < doc->ndiffwords; j++){
			w = doc->indexset[j];
			c = bow[w];
			if (c == 0) continue;
			//printf("Middle : %.2f\n", t.words[w]);
			tw = topic[w];
			if ((tw + c) < lgammalimit)
				loglikelihood += table[tw + c] - table[tw];
			else
				loglikelihood += lgamma(etal + tw + c) - lgamma(etal + tw);
		}
		loglikelihood += gammatable.get(etasum[l],t.nwords);
		loglikelihood -= gammatable.get(etasum[l],doc->levelsum[l] + t.nwords);
		t.loglikelihood = loglikelihood;
	}
}



void hldaTree::newTopic(list<Topic>::iterator t){
	topics.emplace_back(t,nextid); // Empty
	nextid++;
	ntopics++;
}


void hldaTree::ncrp()
{
	for (auto& t : topics)
	{
		if (t.level == 0)
			t.ncrp = 0; // In logarithm 
		else
		{
			t.ncrp = log ((t.ndocs) / (t.parent->ndocs + gamma[t.level - 1]));
		}
	}
}


void hldaTree::cum(vector<Doc>::iterator doc)
{
	list<Topic>::iterator it;
	for (auto& t : topics){
		t.cumlog = t.loglikelihood + t.ncrp;
		it = t.parent;
		for (auto l = t.level - 1; l >= 0; l--)
		{
			t.cumlog += it->loglikelihood + it->ncrp;
			it = it->parent;
		}
		for (auto l = t.level + 1; l < depth; l++)
			t.cumlog += doc->emptyscore[l];

		if (t.level < (depth-1)) // Internal node creates a new branch
			t.cumlog += log(gamma[t.level] / (t.ndocs + gamma[t.level]));
	}
}


list<Topic>::iterator hldaTree::sample()
{
	double maxp=-INFINITY;
	double sump=0;
	double s = urand();

	for (auto& t : topics){
		if (maxp < t.cumlog && t.restriction)
			maxp = t.cumlog;
	}

	for (auto& t : topics){
		t.cumlog = t.cumlog - maxp;
	}

	for (auto& t : topics){
		if (t.restriction)
			t.prob = exp(t.cumlog);
		else
			t.prob = 0;
		sump += t.prob;
	}

	for (auto& t : topics){
		t.prob = t.prob / sump;
	}

	list<Topic>::iterator t;
	for (t = topics.begin();t!=topics.end();t++){
		if (t->prob > s)
			return t;
		s = s - t->prob;
	}

}


void hldaTree::addToPath(vector<Doc>::iterator d, list<Topic>::iterator topic,int heldout)
{
	list<Topic>::iterator nt=topic;
	while (nt->level != depth-1) // Internal node create children
	{
		newTopic(nt);
		nt = --topics.end(); // Last item
	}
	int j,at,bet;
	while (true) {
		nt->ndocs++;
			for ( j = 0; j < d->ndiffwords; j++)
			{
				bet = d->indexset[j];
				at = d->bow.data[nt->level*nterms + bet];
				nt->words[bet] += at ; // Do it more efficiently !!!!
			}
		
		nt->nwords += d->levelsum[nt->level];
		d->path[nt->level] = nt;
		if (nt->level==(depth-1) && heldout == 0)
			nt->docs.push_back(d);

		if (nt->level > 0)
			nt = nt->parent;
		else break;
	}
	d->validpath = 1;

}

void hldaTree::removeFromPath(vector<Doc>::iterator d,int heldout)
{
	int count,w;
	double* dtopic;
	double* bow;
	for (auto i = depth-1; i >= 0; i--)
	{
		d->path[i]->ndocs--;

		if (i == (depth - 1) && heldout == 0) //Leaf node
		{
			auto it = d->path[i]->docs.begin();
			int found = 0;
			if ((*it) == d) 
				d->path[i]->docs.erase(it); // Be careful docs should came in the same order
			else {
				for (; it != d->path[i]->docs.end(); it++) // Search
				if (*it == d) {
					d->path[i]->docs.erase(it);
					found = 1;
					break;
				}
				if (!found) throw 1;   // Error
			}
			
		}


		if (d->path[i]->ndocs == 0 && i!=0)
			removeTopic(d->path[i]);
		else
		{
			bow = d->bow.data+i*nterms;
			dtopic = d->path[i]->words.data;
			// d->path[i]->words = d->path[i]->words - d->bow[i];
			for (int j = 0; j < d->ndiffwords; j++)
			{
				w = d->indexset[j];
				count = bow[w];
				if (count > 0)
					dtopic[w] = dtopic[w] - count; // Do it more efficiently !!!!
			}
			d->path[i]->nwords -= d->levelsum[i];
		}
	}
	d->validpath = 0;
}

void hldaTree::removeTopic(list<Topic>::iterator topic)
{
	topics.erase(topic);
	ntopics--;
}