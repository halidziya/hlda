#include "Sampler.h"
#include "params.h"

Sampler::Sampler() : state()
{
}

Sampler::Sampler(string filename) : state(filename)
{
}

Sampler::Sampler(string datafile, string confile) : state(datafile, confile)
{
}	


Sampler::Sampler(string datafile, string confile, string testfile) : state(datafile, confile, testfile)
{
}

Sampler::~Sampler()
{
}


void Sampler::initialize()
{

	// generator.seed(time(NULL));
	// srand(time(NULL));
	
	auto gammavector = [](double add,int size){
		Vector res(size);
		for (int i = 0; i < size; i++)
			res[i] = lgamma(i + add);
		return res;
	};


	// Precompute gamma values
	lgammalimit = state.docs[0].length * state.ndocs;
	for (auto i = 0; i < depth; i++)
	{
		gammatable[eta[i]] = gammavector(eta[i], lgammalimit);
		gammatable[etasum[i]] = gammavector(etasum[i], lgammalimit);
		gammatable[alpha[i]] = gammavector(alpha[i], lgammalimit);
		if (i < depth - 1)
		{
			gammatable[gamma[i]] = gammavector(gamma[i], lgammalimit);
			gammatable[gamma[i] * 100] = gammavector(gamma[i], lgammalimit);
		}
	}
	gammatable[0] = gammavector(0, lgammalimit);
	gammatable[alpha.sum()] = gammavector(alpha.sum(), lgammalimit);
	
	for (auto i = 0; i < inititer; i++)
	{
		for (Doc& d : state.docs)
			d.samplezfromPrior(); 
	
		
		for (auto d = state.docs.begin(); d != state.docs.end(); d++)
		{
			samplePath(d, 1);
			d->sampleLevels();
			samplePath(d, 0);
			d->sampleLevels();
		}
		cout << "Init: " << i << endl;
	
		state.likelihood(); // Calculate likelihood and store into score

		if (i == 0 || state.score > beststate.score)
		{
			beststate = state;
		}

		for (auto d = state.docs.begin(); d != state.docs.end(); d++)
			state.tree.removeFromPath(d);
		
	}
	state = beststate;

}


double Sampler::samplePath(vector<Doc>::iterator d,State& state,int init,int heldout)
{
	if (!init)
		state.tree.removeFromPath(d,heldout);
	d->emptyScore();
	state.tree.likelihoods(d);
	state.tree.ncrp();
	state.tree.cum(d);
	auto t = state.tree.sample();
	state.tree.addToPath(d, t,heldout);
	return log(t->prob);
}

double Sampler::immitatePath(vector<Doc>::iterator d,vector<Doc>::iterator copy , State& state)
{
	state.tree.removeFromPath(d);
	d->emptyScore();
	state.tree.likelihoods(d);
	state.tree.ncrp();
	state.tree.cum(d);
	// auto t = state.tree.sample();
	// find corresponding node
	
	list<Topic>::iterator t;
	for (int i = 0; i < depth;i++)
	if (copy->path[i]->restriction)
		t = copy->path[i]->copy;
	state.tree.addToPath(d, t);
	for (int i = 0; i < depth; i++)
	{
		copy->path[i]->restriction = 1;
		copy->path[i]->copy = d->path[i];
	}
	return log(t->prob);
}

double Sampler::samplePath(vector<Doc>::iterator d, int init,int heldout)
{
	return samplePath(d, state, init,heldout);
}

void Sampler::gibbs()
{
	for (auto d = state.docs.begin(); d != state.docs.end(); d++)
	{
		samplePath(d, 0);
		d->sampleLevels();
	}
	//cout << "Score : " << score << endl;
}

// Modify model to keep list only in the leaf nodes. 
// It will be easy to remove from those lists
// Collect list in order, then call tree.removepath for documents. Since order satisfied only remove the begin()
// List also be used in sampling. After i and j put restrict others to be either paths or in that path 

// How to identify leaf nodes in restriction , if parent is restricted it is also restricted
// Combine list of restricted ones
void Sampler::mh()
{
	int s, k, r, l;
	list<vector<Doc>::iterator>::iterator di;
	vector<Doc>::iterator mi, mj, si, sj, i, j;
	list<Topic>::iterator st, tit, mtit, stit,otit;
	list<vector<Doc>::iterator> docs,mergedocs,splitdocs,test;
	if (smiter == 0) return;
	l = rand() % (depth - 1); // No SM in last level
	s = 0;
	for (auto& t : state.tree.topics)
	{
		if (t.level == l)
		{
			s += (t.ndocs * (t.ndocs - 1)) / 2; // Possible pairs
		}
	}

	r = rand() % s; // Select node based on pairs, uniformly on pairs
	for (tit = state.tree.topics.begin(); tit != state.tree.topics.end(); tit++)
	{
		if (tit->level == l)
		{
			r = r - (tit->ndocs * (tit->ndocs - 1)) / 2;
			if (r <= 0) break;
		}
	}


	for (auto& t : state.tree.topics)
		t.restriction = 0;
	tit->restriction = 1;
	st = tit;
	// Collect document list in that node
	for (tit = (++state.tree.topics.begin()); tit != state.tree.topics.end(); tit++)
	{
		if (tit->parent->restriction == 1)
		{
			tit->restriction = 1;
		}
		if (tit->level == (depth - 1) && tit->restriction == 1)
			docs.insert(docs.end(), tit->docs.begin(), tit->docs.end());
	}
	//printf("Docs : %d,%d \n",st->ndocs, docs.size());
	
	r = rand() % docs.size();
	for (di = docs.begin(); (di != docs.end()) && (r != 0); di++, r--);
	i = *di;
	docs.erase(di);

	r = rand() % docs.size();
	for (di = docs.begin(); (di != docs.end()) && (r != 0); di++, r--);
	j = *di;

	docs.clear();
	for (auto& t : state.tree.topics)
		t.restriction = 0;
	i->path[l + 1]->restriction = 1;
	j->path[l + 1]->restriction = 1;
	for (tit = (++state.tree.topics.begin()); tit != state.tree.topics.end(); tit++)
	{
		if (tit->parent->restriction == 1)
		{
			tit->restriction = 1;
		}
		if (tit->level == (depth - 1) && tit->restriction == 1)
		{
			for (auto dit = tit->docs.begin(); dit != tit->docs.end(); dit++)
			if ((*dit != i) && (*dit != j)) // i and j excluded
			{
				docs.insert(docs.end(), *dit);
			}
		}
	}

	for (auto& t : state.tree.topics) // Keep the original copy
		t.restriction = 1;
	//printf("Level: %d ,selected pairs: %d %d\nsplit : %d\nndocs:%d %d\n", l, i->id, j->id, i->path[l + 1] == j->path[l + 1], st->ndocs, docs.size());

	// Merge launch
	mergelaunch = state;

	for (auto& t : mergelaunch.tree.topics)
		t.restriction = 0;

	mtit = st->copy; // Selected node
	mi = i->copy;
	mj = j->copy;
	mi->path[l + 1]->restriction = 1; // Docs in child node i
	mj->path[l + 1]->restriction = 1; // Docs in child node j
	for (tit = (++mergelaunch.tree.topics.begin()), otit = (++state.tree.topics.begin()); tit != mergelaunch.tree.topics.end(); tit++,otit++)
	{
		if (tit->parent->restriction == 1)
		{
			tit->restriction = 1;
		}
	}
	for (auto& di : docs)
		mergedocs.insert(mergedocs.end(), di->copy);


	for (auto& a : mergedocs) // Randomly assign level allocations
	{
		mergelaunch.tree.removeFromPath(a);
		//a->samplezfromPrior();
		//a->allLevel2();
	}


	mergelaunch.tree.removeFromPath(mj); // Remove from merged node
	for (k = 0; k < depth; k++)			 // Sample levels of mi in its path
	{
		mi->path[k]->words   = mi->path[k]->words - mi->bow[k];
		mi->path[k]->nwords -= mi->levelsum[k];
	}
	//mi->samplezfromPrior();
	//mi->allLevel2();
	for (k = 0; k < depth; k++)
	{
		mi->path[k]->nwords += mi->levelsum[k];
		mi->path[k]->words  = mi->path[k]->words + mi->bow[k];
	}

	//mj->samplezfromPrior();
	//mj->allLevel2();
	// Merge if not merged
	samplePath(mj,mergelaunch,1);

	mi->sampleLevels();
	mj->sampleLevels();
	for (auto& a : mergedocs) // Randomly assign level allocations
	{
		samplePath(a, mergelaunch,1);
		a->sampleLevels();
	}

	// Merge launch initial is ready
	// Split launch
	
	splitlaunch = state;
	stit = st->copy; // Selected node
	
	for (auto& t : splitlaunch.tree.topics)
		t.restriction = 0;

	si = i->copy;
	sj = j->copy;
	si->path[l + 1]->restriction = 1; // Docs in child node i
	sj->path[l + 1]->restriction = 1; // Docs in child node j
	

	// Mark children nodes too
	for (tit = (++splitlaunch.tree.topics.begin()); tit != splitlaunch.tree.topics.end(); tit++)
	{
		if (tit->parent->restriction == 1)
		{
			tit->restriction = 1;
		}
	}

	for (auto& di : docs)
		splitdocs.insert(splitdocs.end(), di->copy);

	for (auto& a : splitdocs) // Randomly assign level allocations
	{
		splitlaunch.tree.removeFromPath(a); // Corresponding doc in 
		//a->samplezfromPrior();
		//a->allLevel2();
	}


	// Remove words
	for (k = 0; k < depth; k++)
	{
		si->path[k]->words = si->path[k]->words - si->bow[k];
		sj->path[k]->words = sj->path[k]->words - sj->bow[k];
		si->path[k]->nwords -= si->levelsum[k];
		sj->path[k]->nwords -= sj->levelsum[k];
	}
	//si->samplezfromPrior();
	//sj->samplezfromPrior();
	//si->allLevel2();
	//sj->allLevel2();
	// Add words
	for (k = 0; k < depth; k++)
	{
		si->path[k]->words = si->path[k]->words + si->bow[k];
		si->path[k]->nwords += si->levelsum[k];
		sj->path[k]->words = sj->path[k]->words + sj->bow[k];
		sj->path[k]->nwords += sj->levelsum[k];

	}

	// Split if not splitted
	if (si->path[l + 1] == sj->path[l + 1])
	{
		splitlaunch.tree.removeFromPath(sj); // Remove from merged node
		splitlaunch.tree.addToPath(sj,stit); // Add to another branch
	}

	si->sampleLevels();
	sj->sampleLevels();
	for (auto& a : splitdocs) // Randomly assign level allocations
	{
		samplePath(a, splitlaunch, 1);
		a->sampleLevels();
	}

	printf("Sweet with #docs: %d in level %d\n", splitdocs.size(),l);
	// SWEETENING
	for (k = 0; k < smiter; k++)
	{

		if (i->path[l + 1] == j->path[l + 1]){
			si->sampleLevels();
			sj->sampleLevels();

			// Sample in their own subtree
			for (tit = (++splitlaunch.tree.topics.begin()); tit != splitlaunch.tree.topics.end(); tit++)
				tit->restriction = 0;
			si->path[l + 1]->restriction = 1; // Docs in child node i
			for (tit = (++splitlaunch.tree.topics.begin()); tit != splitlaunch.tree.topics.end(); tit++)
			if (tit->parent->restriction == 1)
				tit->restriction = 1;
			if (si->path[l + 1]->ndocs != 1) // Only one in restricted path
				samplePath(si, splitlaunch, 0);

			// Sample in their own subtree
			for (tit = (++splitlaunch.tree.topics.begin()); tit != splitlaunch.tree.topics.end(); tit++)
				tit->restriction = 0;
			sj->path[l + 1]->restriction = 1; // Docs in child node j
			for (tit = (++splitlaunch.tree.topics.begin()); tit != splitlaunch.tree.topics.end(); tit++)
			if (tit->parent->restriction == 1)
				tit->restriction = 1;
			if (sj->path[l + 1]->ndocs != 1) // Only one in restricted path
				samplePath(sj, splitlaunch, 0);


			// Other docs can go either path
			for (tit = (++splitlaunch.tree.topics.begin()); tit != splitlaunch.tree.topics.end(); tit++)
				tit->restriction = 0;
			sj->path[l + 1]->restriction = 1;
			si->path[l + 1]->restriction = 1;
			for (tit = (++splitlaunch.tree.topics.begin()); tit != splitlaunch.tree.topics.end(); tit++)
			if (tit->parent->restriction == 1)
				tit->restriction = 1;

			for (auto& a : splitdocs)
			{
				a->sampleLevels();
				samplePath(a, splitlaunch, 0);
			}
		}
		if (i->path[l + 1] != j->path[l + 1]){
			mi->sampleLevels();
			mj->sampleLevels();
			samplePath(mj, mergelaunch, 0);
			samplePath(mi, mergelaunch, 0);
			for (auto& a : mergedocs)
			{
				a->sampleLevels();
				samplePath(a, mergelaunch, 0);
			}
		}
	}
	
	state.likelihood();
	// Final Scan
	double qup=0, qdown=0;
	if (i->path[l + 1] == j->path[l + 1]) // Split scan
	{
		
		qdown = si->sampleLevels(1);
		qdown += sj->sampleLevels(1);

		for (auto& a : splitdocs)
		{
			qdown += a->sampleLevels(1);
			//qdown += samplePath(a, splitlaunch, 0);
		}

		
		qup  = i->immitateLevels(i);
		qup += j->immitateLevels(j);

		for (tit = state.tree.topics.begin(); tit != state.tree.topics.end(); tit++)
			tit->restriction = 0;
		for (int li = 0; li < depth; li++){
			i->path[li]->copy = mi->path[li];
			i->path[li]->restriction = 1;
		}


		// Self transition , just transition to itself so path is same
		immitatePath(mj, j, mergelaunch);
		for (auto md = mergedocs.begin(), od = docs.begin(); md != mergedocs.end(); md++, od++)
		{
			immitatePath(*md, *od, mergelaunch);
		}
		for (auto md = mergedocs.begin(), od = docs.begin(); md != mergedocs.end(); md++, od++)
		{
			qup += (*od)->immitateLevels(*od);
		}

		for (tit = state.tree.topics.begin(); tit != state.tree.topics.end(); tit++)
			tit->restriction = 1;

		//printf("Split Propose"); // %.2f\n", splitlaunch.likelihood() + qup - qdown - state.score);		
		if (exp(splitlaunch.likelihood() - state.score  ) > urand())
		{
			printf("Split Propose %.2f , State : %.2f , QUP : %.2f , QDOWN : %.2f\n", splitlaunch.likelihood(), state.score, qup, qdown);
			printf("level : %d, ntopics %d->%d\n", l, state.tree.ntopics, splitlaunch.tree.ntopics);
			// printf("O=doc : %d , %d %d %d  doc2: %d , %d %d %d\n", i->id, i->path[0]->ndocs, i->path[1]->ndocs, i->path[2]->ndocs, j->id, j->path[0]->ndocs, j->path[1]->ndocs, j->path[2]->ndocs);
			// printf("M=doc : %d , %d %d %d  doc2: %d , %d %d %d\n", si->id, si->path[0]->ndocs, si->path[1]->ndocs, si->path[2]->ndocs, sj->id, sj->path[0]->ndocs, sj->path[1]->ndocs, sj->path[2]->ndocs);
			for (auto& t : splitlaunch.tree.topics)
				t.restriction = 1;
			printf("=============\nSPLIT ACCEPTED\n=============\n%d %d\n%.2f %.2f\n", state.tree.ntopics, splitlaunch.tree.ntopics,state.score,splitlaunch.score);
			state = splitlaunch;
			//system("pause");
		}
	}
	else
	{
		qdown = mi->sampleLevels(1);
		qdown += mj->sampleLevels(1);
		//qdown += samplePath(mj, mergelaunch, 0);
		for (auto& a : mergedocs)
		{
			qdown += a->sampleLevels(1);
			//qdown += samplePath(a, mergelaunch, 0);
		}

		for (tit = state.tree.topics.begin(); tit != state.tree.topics.end(); tit++)
			tit->restriction = 0;
		for (int li = 0; li < depth; li++){
			i->path[li]->copy = si->path[li];
			j->path[li]->copy = sj->path[li];
			i->path[li]->restriction = 1;
			j->path[li]->restriction = 1;
		}


		for (auto md = splitdocs.begin(), od = docs.begin(); md != splitdocs.end(); md++, od++)
		{
			immitatePath(*md, *od, splitlaunch);
		}

		qup = i->immitateLevels(i);
		qup += j->immitateLevels(j);
		for (auto md = splitdocs.begin(), od = docs.begin(); md != splitdocs.end(); md++, od++)
			qup += (*od)->immitateLevels(*od);

		for (tit = state.tree.topics.begin(); tit != state.tree.topics.end(); tit++)
			tit->restriction = 1;


		
		
		if (exp(mergelaunch.likelihood() - state.score  ) > urand()) // + qup - qdown
		{
			printf("Merge Propose %.2f , State : %.2f , QUP : %.2f , QDOWN : %.2f\n", mergelaunch.likelihood(), state.score, qup, qdown);
			printf("level : %d, ntopics %d->%d\n", l, state.tree.ntopics, mergelaunch.tree.ntopics);
			// printf("O=doc : %d , %d %d %d  doc2: %d , %d %d %d\n", i->id, i->path[0]->ndocs, i->path[1]->ndocs, i->path[2]->ndocs, j->id, j->path[0]->ndocs, j->path[1]->ndocs, j->path[2]->ndocs);
			// printf("M=doc : %d , %d %d %d  doc2: %d , %d %d %d\n", mi->id, mi->path[0]->ndocs, mi->path[1]->ndocs, mi->path[2]->ndocs, mj->id, mj->path[0]->ndocs, mj->path[1]->ndocs, mj->path[2]->ndocs);
			for (auto& t : mergelaunch.tree.topics)
				t.restriction = 1;

			printf("=============\nMERGE ACCEPTED\n=============\n%d %d\n%.2f %.2f\n", state.tree.ntopics, mergelaunch.tree.ntopics, state.score, mergelaunch.score);
			state = mergelaunch;
			
		}
		//system("pause");
		
	}


	// Proposal Calculations


	// Likelihoods
	//printf("Original Score : %.2f %d\n", state.score,state.tree.ntopics);
	/*printf("Merge Score : %.2f\n", mergelaunch.likelihood());
	printf("Split Score : %.2f\n", splitlaunch.likelihood());*/

	
	// Accept or Reject
}



double Sampler::heldoutScore() // Based on Blei's implementation
{
	int burnin = 200;
	int lag = 10;
	int hiter = 500;
	double trainscore = state.likelihood();
	double oldeta = state.etaScore();
	double score = 0,s=0;
	int count = 0;
	int nwords = 0;
	double alphasum = alpha.sum();

/*	for (auto d = state.docs.begin(); d != state.docs.end(); d++)
	{
		d->sampleunifheld();
		//d->allLevel2();
		samplePath(d, 0, 1);
		d->sampleLevelsheld();
		nwords += d->nheld;
		for (int j = 0; j < d->nheld; j++)
			d->like[j] = 0;
	}

	for (auto i = 0; i <= hiter; i++)
	{
		//Iterate gibbs state
		for (auto d = state.docs.begin(); d != state.docs.end(); d++)
		{
			samplePath(d, 0, 1);
			d->sampleLevelsheld();
		}
		if ((i > burnin) && ((i % lag) == 0))
		{
			for (auto d = state.docs.begin(); d != state.docs.end(); d++)
			{
				for (int j = 0; j < d->nheld; j++)
				{
					int w = d->heldout[j];
					s = 0;
					for (int l = 0; l < depth; l++)
					{
						s += ((d->levelsum[l] + alpha[l]) / (d->length + alphasum)) * ((d->path[l]->words[w] + eta[l]) / (d->path[l]->nwords + etasum[l]));
					}
					d->like[j] += s;
				}
			}
		}
	}

	for (auto d = state.docs.begin(); d != state.docs.end(); d++)
	{
		for (int j = 0; j < d->nheld; j++)
		{
			d->like[j] /= ((hiter - burnin) / lag);
			score += log(d->like[j]);
		}
	}*/
	
	// cout << "Likelihood " << state.likelihood() << "Eta: " << state.etaScore() << "Gamma: " << state.crpScore() << "Alpha: " << state.alphaScore() << endl;
	if (nheldout > 0)
	{
		// Initialize
		for (auto d = state.heldoutdocs.begin(); d != state.heldoutdocs.end(); d++)
		{
			d->sampleunif();
			//d->allLevel2();
			samplePath(d, 1,1);
			d->sampleLevels();
			nwords += d->length;
			d->like = -INFINITY;
		}

		score = 0;
		for (auto i = 0; i < hiter; i++)
		{
			//Iterate gibbs state
			for (auto d = state.heldoutdocs.begin(); d != state.heldoutdocs.end(); d++)
			{
				samplePath(d, 0,1);
				d->sampleLevels();
			}
			if ((i >= burnin) && ((i % lag) == 0))
			{
				count++;
				for (auto d = state.heldoutdocs.begin(); d != state.heldoutdocs.end(); d++)
				{
					
					score = 0;
					for (int j = 0; j < d->length;j++)
					{
						int w = d->words[j];
						s = 0;
						for (int l = 0; l < depth; l++)
						{
						///int l = d->z[j];
						s += ((d->levelsum[l] + alpha[l]) / (d->length + alphasum)) * ((d->path[l]->words[w] + eta[l]) / (d->path[l]->nwords + etasum[l]));
						}
						score += log(s);
					}
					
					// Logsum
					if (d->like > score)
						d->like = d->like + log(1 + exp(score - d->like));
					else
						d->like = score + log(1 + exp(d->like - score));
					
				}

			}

		}
		
		score = 0;
		for (auto d = state.heldoutdocs.begin(); d != state.heldoutdocs.end(); d++)
		{
				d->like -= log(count);
				score += (d->like);
		}
	}
	return score/(nwords);
}


void Sampler::parallelmh(){
	if (smiter == 0) return;
	int nthreads =  thread::hardware_concurrency();
	Vector acceptance(nthreads);
	Vector splitormerge(nthreads);
	Vector qups(nthreads);
	Vector qdowns(nthreads);
	Vector ls(nthreads);
	Vector rs(nthreads);
	Vector dis(nthreads);
	Vector djs(nthreads);
	for (int th_i = 0; th_i < nthreads; th_i++)
	{
		copies[th_i] = state;
		// Different random numbers for each thread
		ls[th_i] = rand() % (depth - 1);
		rs[th_i] = rand();
		dis[th_i] = rand();
		djs[th_i] = rand();
	}
	state.likelihood();
	parallel_for(size_t(0), (size_t)nthreads, [&](size_t th_i)
	{
		buffer.threadid = th_i;
		absbuffer.threadid = th_i;
		State& copy = copies[th_i];
		int l = ls[th_i];
		int s = 0, r = 0;
		int k;
		srand(th_i+time(NULL));
		seedMT(th_i + time(NULL));
		copy.likelihood();
		list<Topic>::iterator st, tit, mtit, stit, otit;
		vector<Doc>::iterator mi, mj, si, sj, i, j;
		list<vector<Doc>::iterator>::iterator di;
		list<vector<Doc>::iterator> docs, mergedocs, splitdocs, test;
		for (auto& t : copy.tree.topics)
		{
			if (t.level == l)
			{
				s += (t.ndocs * (t.ndocs - 1)) / 2; // Possible pairs
			}
		}
		r = (int)rs[th_i] % s; // Select node based on pairs, uniformly on pairs
		for (tit = copy.tree.topics.begin(); tit != copy.tree.topics.end(); tit++)
		{
			if (tit->level == l)
			{
				r = r - (tit->ndocs * (tit->ndocs - 1)) / 2;
				if (r <= 0) break;
			}
		}

		for (auto& t : copy.tree.topics)
			t.restriction = 0;
		tit->restriction = 1;
		st = tit;


		for (tit = (++copy.tree.topics.begin()); tit != copy.tree.topics.end(); tit++)
		{
			if (tit->parent->restriction == 1)
			{
				tit->restriction = 1;
			}
			if (tit->level == (depth - 1) && tit->restriction == 1)
				docs.insert(docs.end(), tit->docs.begin(), tit->docs.end());
		}

		r = (int)dis[th_i] % docs.size();
		for (di = docs.begin(); (di != docs.end()) && (r != 0); di++, r--);
		i = *di;
		docs.erase(di);

		r = (int)djs[th_i] % docs.size();
		for (di = docs.begin(); (di != docs.end()) && (r != 0); di++, r--);
		j = *di;

		docs.clear();
		for (auto& t : copy.tree.topics)
			t.restriction = 0;
		i->path[l + 1]->restriction = 1;
		j->path[l + 1]->restriction = 1;

		for (tit = (++copy.tree.topics.begin()); tit != copy.tree.topics.end(); tit++)
		{
			if (tit->parent->restriction == 1)
			{
				tit->restriction = 1;
			}
			if (tit->level == (depth - 1) && tit->restriction == 1)
			{
				for (auto dit = tit->docs.begin(); dit != tit->docs.end(); dit++)
				if ((*dit != i) && (*dit != j)) // i and j excluded
				{
					docs.insert(docs.end(), *dit);
				}
			}
		}
		
		for (auto& t : copy.tree.topics) // Keep the original copy
			t.restriction = 1;

		mergelaunches[th_i] = copy;
		for (auto& t : mergelaunches[th_i].tree.topics)
			t.restriction = 0;


		mtit = st->copy; // Selected node
		mi = i->copy;
		mj = j->copy;
		mi->path[l + 1]->restriction = 1; // Docs in child node i
		mj->path[l + 1]->restriction = 1; // Docs in child node j

		for (tit = (++mergelaunches[th_i].tree.topics.begin()), otit = (++copy.tree.topics.begin()); tit != mergelaunches[th_i].tree.topics.end(); tit++, otit++)
		{
			if (tit->parent->restriction == 1)
			{
				tit->restriction = 1;
			}
		}


		for (auto& di : docs)
			mergedocs.insert(mergedocs.end(), di->copy);

		for (auto& a : mergedocs) // Randomly assign level allocations
		{
			mergelaunches[th_i].tree.removeFromPath(a);
			a->samplezfromPrior();
		}


		mergelaunches[th_i].tree.removeFromPath(mj); // Remove from merged node

		for (k = 0; k < depth; k++)			 // Sample levels of mi in its path
		{
			// mi->path[k]->words = mi->path[k]->words - mi->bow[k];
			//for (int& w : mi->indexset)
			//	mi->path[k]->words.data[w] -= mi->bow.data[k*nterms + w];
			// mi->path[k]->words = mi->path[k]->words - mi->bow[k];
			mi->path[k]->words = mi->path[k]->words - mi->bow[k];
			mi->path[k]->nwords -= mi->levelsum[k];
		}
		mi->samplezfromPrior();
		for (k = 0; k < depth; k++)
		{
			//for (int& w : mi->indexset)
			//	mi->path[k]->words.data[w] += mi->bow.data[k*nterms + w];
			mi->path[k]->words = mi->path[k]->words + mi->bow[k];
			mi->path[k]->nwords += mi->levelsum[k];
		}

		mj->samplezfromPrior();
		// Merge if not merged
		samplePath(mj, mergelaunches[th_i], 1);
		mi->sampleLevels();
		mj->sampleLevels();
		for (auto& a : mergedocs) // Randomly assign level allocations
		{
			samplePath(a, mergelaunches[th_i], 1);
			a->sampleLevels();
		}

		// Merge launch initial is ready
		// Split launch

		splitlaunches[th_i] = copy;
		stit = st->copy; // Selected node

		for (auto& t : splitlaunches[th_i].tree.topics)
			t.restriction = 0;

		si = i->copy;
		sj = j->copy;
		si->path[l + 1]->restriction = 1; // Docs in child node i
		sj->path[l + 1]->restriction = 1; // Docs in child node j

		// Mark children nodes too
		for (tit = (++splitlaunches[th_i].tree.topics.begin()); tit != splitlaunches[th_i].tree.topics.end(); tit++)
		{
			if (tit->parent->restriction == 1)
			{
				tit->restriction = 1;
			}
		}

		for (auto& di : docs)
			splitdocs.insert(splitdocs.end(), di->copy);

		for (auto& a : splitdocs) // Randomly assign level allocations
		{
			splitlaunches[th_i].tree.removeFromPath(a); // Corresponding doc in 
			a->samplezfromPrior();
		}

		// Remove words
		for (k = 0; k < depth; k++)
		{
		/*	for (int& w : si->indexset)
				si->path[k]->words.data[w] -= si->bow.data[k*nterms + w];
			for (int& w : sj->indexset)
				sj->path[k]->words.data[w] -= sj->bow.data[k*nterms + w];*/
			si->path[k]->words = si->path[k]->words - si->bow[k];
			sj->path[k]->words = sj->path[k]->words - sj->bow[k];
			si->path[k]->nwords -= si->levelsum[k];
			sj->path[k]->nwords -= sj->levelsum[k];
		}

		si->samplezfromPrior();
		sj->samplezfromPrior();

		// Add words
		for (k = 0; k < depth; k++)
		{
			/*for (int& w : si->indexset)
				si->path[k]->words.data[w] += si->bow.data[k*nterms + w];
			for (int& w : sj->indexset)
				sj->path[k]->words.data[w] += sj->bow.data[k*nterms + w];*/
			si->path[k]->words = si->path[k]->words + si->bow[k];
			si->path[k]->nwords += si->levelsum[k];
			sj->path[k]->words = sj->path[k]->words + sj->bow[k];
			sj->path[k]->nwords += sj->levelsum[k];
		}

		// Split if not splitted
		if (si->path[l + 1] == sj->path[l + 1])
		{
			splitlaunches[th_i].tree.removeFromPath(sj); // Remove from merged node
			splitlaunches[th_i].tree.addToPath(sj, stit); // Add to another branch
		}

		si->sampleLevels();
		sj->sampleLevels();
		for (auto& a : splitdocs) // Randomly assign level allocations
		{
			samplePath(a, splitlaunches[th_i], 1);
			a->sampleLevels();
		}

		
		for (k = 0; k < smiter; k++)
		{
			if (i->path[l + 1] == j->path[l + 1]){

				si->sampleLevels();
				sj->sampleLevels();

				// Sample in their own subtree
				for (tit = (++splitlaunches[th_i].tree.topics.begin()); tit != splitlaunches[th_i].tree.topics.end(); tit++)
					tit->restriction = 0;
				si->path[l + 1]->restriction = 1; // Docs in child node i
				for (tit = (++splitlaunches[th_i].tree.topics.begin()); tit != splitlaunches[th_i].tree.topics.end(); tit++)
				if (tit->parent->restriction == 1)
					tit->restriction = 1;
				if (si->path[l + 1]->ndocs != 1) // Only one in restricted path
					samplePath(si, splitlaunches[th_i], 0);

				// Sample in their own subtree
				for (tit = (++splitlaunches[th_i].tree.topics.begin()); tit != splitlaunches[th_i].tree.topics.end(); tit++)
					tit->restriction = 0;
				sj->path[l + 1]->restriction = 1; // Docs in child node j
				for (tit = (++splitlaunches[th_i].tree.topics.begin()); tit != splitlaunches[th_i].tree.topics.end(); tit++)
				if (tit->parent->restriction == 1)
					tit->restriction = 1;
				if (sj->path[l + 1]->ndocs != 1) // Only one in restricted path
					samplePath(sj, splitlaunches[th_i], 0);


				// Other docs can go either path
				for (tit = (++splitlaunches[th_i].tree.topics.begin()); tit != splitlaunches[th_i].tree.topics.end(); tit++)
					tit->restriction = 0;
				sj->path[l + 1]->restriction = 1;
				si->path[l + 1]->restriction = 1;
				for (tit = (++splitlaunches[th_i].tree.topics.begin()); tit != splitlaunches[th_i].tree.topics.end(); tit++)
				if (tit->parent->restriction == 1)
					tit->restriction = 1;

				for (auto& a : splitdocs)
				{
					a->sampleLevels();
					samplePath(a, splitlaunches[th_i], 0);
				}
			}
			if (i->path[l + 1] != j->path[l + 1]){
				mi->sampleLevels();
				mj->sampleLevels();

				samplePath(mj, mergelaunches[th_i], 0);
				samplePath(mi, mergelaunches[th_i], 0);
				for (auto& a : mergedocs)
				{
					a->sampleLevels();
					samplePath(a, mergelaunches[th_i], 0);
				}
			}
		}
	//	if (i->path[l + 1] == j->path[l + 1])
	//	printf("Split sweet with #docs: %d in level %d ntopics: %d\n", splitdocs.size(), l, splitlaunches[th_i].tree.ntopics);
	//	else
	//		printf("Merge sweet with #docs: %d in level %d ntopics: %d\n", splitdocs.size(), l, mergelaunches[th_i].tree.ntopics);
		double qup = 0, qdown = 0;
		if (i->path[l + 1] == j->path[l + 1]) // Split scan
		{
			qdown = si->sampleLevels(1);
			qdown += sj->sampleLevels(1);

			for (auto& a : splitdocs)
			{
				qdown += a->sampleLevels(1);
				//qdown += samplePath(a, splitlaunch, 0);
			}

			qup = i->immitateLevels(i);
			qup += j->immitateLevels(j);

			for (tit = copy.tree.topics.begin(); tit != copy.tree.topics.end(); tit++)
				tit->restriction = 0;
			for (int li = 0; li < depth; li++){
				i->path[li]->copy = mi->path[li];
				i->path[li]->restriction = 1;
			}

			immitatePath(mj, j, mergelaunches[th_i]);
			for (auto md = mergedocs.begin(), od = docs.begin(); md != mergedocs.end(); md++, od++)
			{
				immitatePath(*md, *od, mergelaunches[th_i]);
			}
			for (auto md = mergedocs.begin(), od = docs.begin(); md != mergedocs.end(); md++, od++)
			{
				qup += (*od)->immitateLevels(*od);
			}

			for (tit = copy.tree.topics.begin(); tit != copy.tree.topics.end(); tit++)
				tit->restriction = 1;
			if ((splitlaunches[th_i].likelihood() - copy.score + qup - qdown) > 0)
				acceptance[th_i] = 1;
			else
				acceptance[th_i] = exp(splitlaunches[th_i].likelihood() - copy.score + qup - qdown);
			splitormerge[th_i] = 0;
			qups[th_i] = qup;
			qdowns[th_i] = qdown;
		}
		else
		{
			qdown = mi->sampleLevels(1);
			qdown += mj->sampleLevels(1);
			for (auto& a : mergedocs)
			{
				qdown += a->sampleLevels(1);
				//qdown += samplePath(a, mergelaunch, 0);
			}
			for (tit = copy.tree.topics.begin(); tit != copy.tree.topics.end(); tit++)
				tit->restriction = 0;
			for (int li = 0; li < depth; li++){
				i->path[li]->copy = si->path[li];
				j->path[li]->copy = sj->path[li];
				i->path[li]->restriction = 1;
				j->path[li]->restriction = 1;
			}
			for (auto md = splitdocs.begin(), od = docs.begin(); md != splitdocs.end(); md++, od++)
			{
				immitatePath(*md, *od, splitlaunches[th_i]);
			}

			qup = i->immitateLevels(i);
			qup += j->immitateLevels(j);
			for (auto md = splitdocs.begin(), od = docs.begin(); md != splitdocs.end(); md++, od++)
				qup += (*od)->immitateLevels(*od);

			if ((mergelaunches[th_i].likelihood() - copy.score + qup - qdown) > 0)
				acceptance[th_i] = 1;
			else
			acceptance[th_i] = exp(mergelaunches[th_i].likelihood() - copy.score + qup - qdown);
			splitormerge[th_i] = 1;
			qups[th_i] = qup;
			qdowns[th_i] = qdown;
		}
	});
	for (int th_i = 0; th_i < nthreads; th_i++)
	{
	//	printf("\n%d. acceptance  %4.2f qdown: %.2f qup: %.2f\n", th_i, acceptance[th_i], qdowns[th_i], qups[th_i]);
	//	if (splitormerge[th_i] == 1)
	//		printf("MERGE %d %d\n%.2f %.2f\n", state.tree.ntopics, mergelaunches[th_i].tree.ntopics, state.score, mergelaunches[th_i].score);
	//	else
	//		printf("SPLIT %d %d\n%.2f %.2f\n", state.tree.ntopics, splitlaunches[th_i].tree.ntopics, state.score, splitlaunches[th_i].score);
		

		if (acceptance[th_i] >urand())
		{
			if (splitormerge[th_i] == 1)
			{
				printf("=============\nMERGE ACCEPTED\n=============\n%d %d\n%.2f %.2f\n", state.tree.ntopics, mergelaunches[th_i].tree.ntopics, state.score, mergelaunches[th_i].score);
				state = mergelaunches[th_i];
			}
			else
			{
				printf("=============\nSPLIT ACCEPTED\n=============\n%d %d\n%.2f %.2f\n", state.tree.ntopics, splitlaunches[th_i].tree.ntopics, state.score, splitlaunches[th_i].score);
				state = splitlaunches[th_i];
			}
			for (auto& t : state.tree.topics)
				t.restriction = 1;
			break;
		}
	}
}