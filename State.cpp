#include "State.h"
#include <algorithm>

State::State()
{
	ndocs = 0;
}

void State::readDocs(string filename,vector<Doc>& docs)
{
	ifstream input(filename);
	string line;
	while (getline(input, line)) // Read docs
	{
		docs.push_back(Doc(line, docs.size()));
	}
	input.close();
}


void State::readConf(string confile)
{
	string s;
	double d;
	int i;

	ifstream input(confile);

	auto readint = [](ifstream& input){
		string s;
		int i;
		input >> s >> i;
		cout << s << "\t" << i << endl;
		return i;
	};

	auto readline = [](ifstream& input,int depth) {
		string s;
		double d;
		int i;
		Vector val;

		input >> s;
		cout << s << "\t";
		val.resize(depth);

		for (i = 0; i < depth; i++)
		{
			input >> d;
			val[i] = d;
			cout << d << " ";
		}

		cout << endl;
		return val;
	};
	
	depth = readint(input);
	eta = readline(input, depth);
	gamma = readline(input, depth - 1);
	alpha = readline(input, depth);
	nterms = readint(input)+1; // base 0 index ,base 1 index both supported
	init_buffer(8, nterms);
	niter = readint(input);
	inititer = readint(input);
	nheldout = readint(input);
	etasum = eta * nterms;
	cout << "EtaSum: ";
	for (i = 0; i < depth; i++)
	{
		cout << etasum[i] << " ";
	}
	cout << endl;
}

State::State(string filename) : State()
{
	readDocs(filename,docs);
	ndocs = docs.size();
}


State::State(string datafile, string confile)
{
	readConf(confile);
	ndocs = 0;
	tree = hldaTree();
	readDocs(datafile,docs);
	heldoutdocs = vector<Doc>(docs.begin(),docs.begin()+nheldout);
	docs.erase(docs.begin(), docs.begin() + nheldout);
	ndocs = docs.size();
}

State::State(string datafile, string confile,string testfile)
{
	readConf(confile);
	ndocs = 0;
	tree = hldaTree();
	readDocs(testfile,heldoutdocs);
	nheldout = heldoutdocs.size();
	readDocs(datafile,docs);
	ndocs = docs.size();
}


State::~State()
{
}



void State::operator=(State& s){
	vector<Doc>::iterator it;
	list<vector<Doc>::iterator>::iterator lit;
	list<Topic>::iterator tit, otit;
	int i, j;
	docs = s.docs;
	heldoutdocs = s.heldoutdocs;
	tree = s.tree;
	ndocs = s.ndocs;
	score = s.score;
	for (i = 0,it=docs.begin(); i < ndocs; i++,it++)
		s.docs[i].copy = it;
	for (tit = tree.topics.begin(), otit = s.tree.topics.begin(); tit != tree.topics.end(); tit++, otit++)
	{
	otit->copy = tit;
	if (tit != tree.topics.begin()) tit->parent = otit->parent->copy;  // Except for root
	for (lit = tit->docs.begin();lit!= tit->docs.end(); lit++) // Copy doc list of topics
	{
		(*lit) = (*lit)->copy;
	}
	}

	for (i = 0, it = docs.begin(); i < ndocs; i++, it++)
		if (it->validpath)
			for (j = 0; j < depth; j++)
				it->path[j] = it->path[j]->copy;

}




double State::crpScore() // nCRP
{
	int i, j;
	double logsum=0;
	for (auto &t : tree.topics)
	{
		if (t.level < (depth - 1)){
			logsum += gammatable[gamma[t.level]][0] - gammatable[gamma[t.level]][t.ndocs];
		}
		if (t.level) logsum += log(gamma[t.parent->level]) + gammatable[0][t.ndocs];
	}
	return logsum;
}

double State::etaScore() // Dirichlet-Multinomial
{
	int j, w, l;
	double etal;
	double logsum = 0,lgetal;
	for (auto &t : tree.topics){
		l = t.level;
		etal = eta[l];
		lgetal = gammatable.get(etal,0);
		for (w = 0; w < nterms; w++){
			if (t.words[w] == 0) continue;
			logsum  += gammatable.get(etal,t.words[w]);
			logsum -= lgetal;
		}
		logsum += gammatable.get(etasum[l],0);
		logsum -= gammatable.get(etasum[l],t.nwords);
		
	}
	
	return logsum;
}


double State::alphaScore() //Dirichlet-Multinomial
{
	int i, j;
	double logsum=0;
	double alphasumlgamma, sumlgammaalpha,alphasum;

	alphasum = alpha.sum();
	alphasumlgamma = gammatable.get(alphasum,0);
	sumlgammaalpha = 0;
	for (j = 0; j<depth; j++)
	{
		sumlgammaalpha += gammatable.get(alpha[j],0);
	}


	logsum = 0;
	for (auto& d : docs)
	{
		logsum += alphasumlgamma - sumlgammaalpha - gammatable.get(alphasum,d.length);
		for (j = 0; j<depth; j++){
			logsum += gammatable.get(alpha[j],d.levelsum[j]);
		}
	}
	return logsum;
}


double State::heldlikelihood()
{
	int i, j;
	double logsum = 0;
	double alphasumlgamma, sumlgammaalpha, alphasum;

	alphasum = alpha.sum();
	alphasumlgamma = gammatable.get(alphasum,0);
	sumlgammaalpha = 0;
	for (j = 0; j<depth; j++)
	{
		sumlgammaalpha += gammatable[alpha[j]][0];
	}


	logsum = 0;
	
	for (auto& d : docs)
	{
		logsum += alphasumlgamma - sumlgammaalpha - gammatable.get(alphasum,d.length);
		for (j = 0; j<depth; j++){
			logsum += gammatable.get(alpha[j],d.levelsum[j]);
		}
	}
	for (auto& d : heldoutdocs)
	{
		logsum += alphasumlgamma - sumlgammaalpha - gammatable.get(alphasum,d.length);
		for (j = 0; j<depth; j++){
			logsum += gammatable.get(alpha[j],d.levelsum[j]);
		}
	}

	logsum += etaScore() + crpScore();
	return logsum;
}
double State::likelihood()
{
	score = etaScore() + crpScore() + alphaScore();
	return score;
}


void State::write(string filename)
{
	ofstream ftree(filename+".tree");
	ofstream ftopics(filename + ".topics");
	ofstream fpath(filename + ".path");
	// ofstream fdocs(filename + ".docs");

	int i = 0;
	for (auto& t : tree.topics)
	{
		t.id = i+1; // Start from 1
		ftree << i << "\t";
		i++;
	}
	ftree << endl;
	for (auto& t : tree.topics){
		ftree << t.level << "\t" ;
	}

	ftree << endl;
	for (auto& t : tree.topics){
		ftree << t.ndocs << "\t";
	}



	ftree << endl << 0 << "\t";
	for (auto& t : tree.topics){
		if (t.level > 0) ftree << t.parent->id << "\t";
	}

	ftree << endl;
	for (auto& t : tree.topics){
		ftree << t.nwords << "\t";
	}
	ftree << endl;


	
	for (auto& t : tree.topics){
		ftree << t.ncrp << "\t";
	}
	ftree << endl;


	for (auto& t : tree.topics){
		for (i = 0; i < nterms; i++)
		{
			ftopics << t.words[i] << "\t";
		}
		ftopics << endl;
	}

	for (auto& d : docs)
	{
		for (i = 0; i < depth; i++)
			fpath << d.path[i]->id << " ";
		fpath << endl;
	}

	ftree.close();
	ftopics.close();
	fpath.close();
	// fdocs.close();
}