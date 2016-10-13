#include <iostream>
#include "Sampler.h"
#include <Windows.h>


int main(int argc, char** args)
{
	debugOn = 1;
	if (argc < 2)
	{
		cout << "Usage: \nhlda trainfile [configuration] [testfile] [nsplitmerge]" << endl;
		system("pause");
		return -1;
	}

	int outputid;
	ifstream lascount("outputs/lastcount.txt");
	lascount >> outputid;
	lascount.close();
	ofstream olastcount("outputs/lastcount.txt");
	olastcount << ++outputid;
	olastcount.close();
	string foldername = ("outputs/out") + to_string(outputid);
	CreateDirectory( foldername.c_str() , NULL);
	ofstream gibbs(foldername + "/gibbs.log");
	


	seedMT(time(NULL));
	srand(time(NULL));
	Sampler s;
	step();
	smiter = 20;
	if (argc > 4)
	{
		s = Sampler(args[1], args[2], args[3]);
		if (argc > 4)
			smiter = atoi(args[4]);
	}
	else if (argc > 3)
	{
		s = Sampler(args[1],args[2]);
		if (argc > 3)
			smiter = atoi(args[3]);
	}
	else
		s = Sampler(args[1]);
		
	printf("SMITER : %d\n", smiter);
	s.initialize();
	s.state.write(foldername + "/first");
	int nwords = 0;
	for (int i = 0; i < s.state.ndocs; i++)
		nwords += s.state.docs[i].length;
	for (int i = 0; i < niter; i++)
	{
		if (i % 20 == 0)cout << "Iter: " << i << " " << s.state.score  << " ntopics: " << s.state.tree.ntopics << endl;
		//printf("ntopic before gibbs: %d\n", s.state.tree.ntopics);
		s.gibbs();
		s.state.likelihood();
		if (s.state.score > s.beststate.score)
			s.beststate = s.state;
		if (i > 2000 )
		{
			//if (i%8==0)
			//s.parallelmh();
			s.mh();
		} 

		s.state.likelihood();
		gibbs << s.state.score << endl;
		if (s.state.score > s.beststate.score)
			s.beststate = s.state;
	}
	cout << "best :" << s.beststate.score;
	s.state.write(foldername + "/last");
	s.beststate.write(foldername + "/best");
	s.state = s.beststate;
	ofstream ids(foldername + "/ids.txt");
	ofstream levels(foldername + "/levels.txt");
	for (auto& d : s.beststate.docs)
	{
		ids << d.id << endl;
		for (int l = 0; l < depth; l++)
			levels << d.levelsum[l] << " ";
		levels << endl;
	}
	ids.close();
	levels.close();
	step();
	gibbs.close();
	string settingname = foldername + "/settings.txt";
	CopyFile(args[2], settingname.c_str(), NULL);
	if (nheldout > 0){
		double score = s.heldoutScore();
		ofstream heldout(foldername + "/heldout." + to_string(smiter) + "." + to_string(score));
		heldout << smiter << endl << score << endl;
		cout << "Held Out Score: " << score << endl;
		heldout.close();
	}
	//system("pause");
	
}