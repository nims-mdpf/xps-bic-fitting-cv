#include "genetic_algorithm.h"

DataObj d;

bool readConfig(InfoGenetic& info){
	// Default setting value
	info.MaxOfLoop = 20;
	info.generation = 400;
	info.individual = 30;
	info.num_island = 60;
	info.gw_alpha = 40.0; // for gaussian std
	info.gw_beta = 0.1;
	info.lw_alpha = 60.0; // for lorentzian HWHM
	info.lw_beta = 0.03;
	info.prob_cross = 0.9;
	info.prob_exchange = 0.3;
	info.prob_mutation = 0.2;
	info.beta_log_upper = 1.25;
	info.beta_log_lower = -1;
	info.bgh_alpha = 2.0;
	info.bgh_beta = 100.0;
	info.bgl_alpha = 2.0;
	info.bgl_beta = 100.0;
	// info.bg_stddev = 6;     // std of bg for both side
	// info.shift_stddev = 0.2;   // not used
	info.intensity_alpha = 1;   // for height 
	info.intensity_beta = 1100;   // for height 
	info.intensity_max = 1e300;
	info.noiseType = "Poisson";
	

    const string fname = "config.txt";
    ifstream ifs(fname);
    string buf, str;
    vector<string> v;
    if(ifs.fail()){
		cout << "Failed to read " << fname << endl;
		return true; 
	}
    while(getline(ifs, buf)){
        //replace(buf.begin(), buf.end(), ':', ' ');
		if (buf.empty() || buf[0] == '#') continue;
		buf.erase(remove(buf.begin(), buf.end(),' '), buf.end());
        stringstream line{buf};
        while(getline(line, str, ':')){ v.push_back(str); }
		if(v[0]=="MaxOfLoop"){ info.MaxOfLoop = stoi(v[1]); }
        if(v[0]=="generation"){ info.generation = stoi(v[1]); }
        if(v[0]=="individual"){ info.individual = stoi(v[1]); }
		if(v[0]=="num_island"){ info.num_island = stoi(v[1]); }
		if(v[0]=="beta_log_upper"){ info.beta_log_upper = stoi(v[1]); }
		if(v[0]=="beta_log_lower"){ info.beta_log_lower = stoi(v[1]); }
        if(v[0]=="gw_alpha"){ info.gw_alpha = stod(v[1]); }
        if(v[0]=="gw_beta"){ info.gw_beta = stod(v[1]); }
        if(v[0]=="lw_alpha"){ info.lw_alpha = stod(v[1]); }
        if(v[0]=="lw_beta"){ info.lw_beta = stod(v[1]); }
        if(v[0]=="bgh_alpha"){ info.bgh_alpha = stod(v[1]); }
        if(v[0]=="bgh_beta"){ info.bgh_beta = stod(v[1]); }
        if(v[0]=="bgl_alpha"){ info.bgl_alpha = stod(v[1]); }
        if(v[0]=="bgl_beta"){ info.bgl_beta = stod(v[1]); }
		// if(v[0]=="bg_stddev"){ info.bg_stddev = stod(v[1]); }
		// if(v[0]=="shift_stddev"){ info.shift_stddev = stod(v[1]); }
		if(v[0]=="intensity_alpha"){ info.intensity_alpha = stod(v[1]); }
		if(v[0]=="intensity_beta"){ info.intensity_beta = stod(v[1]); }
		if(v[0]=="intensity_max"){ info.intensity_max = stod(v[1]); }
		if(v[0]=="prob_exchange"){ info.prob_exchange = stod(v[1]); }
		if(v[0]=="prob_cross"){ info.prob_cross = stod(v[1]); }
		if(v[0]=="prob_mutation"){ info.prob_mutation = stod(v[1]); }
		if(v[0]=="noiseType"){ info.noiseType = v[1]; }
        v.clear();
    }

	cout <<"info.MaxOfLoop       = " << info.MaxOfLoop << endl;
	cout <<"info.generation      = " << info.generation << endl;
	cout <<"info.individual      = " << info.individual << endl;
	cout <<"info.num_island      = " << info.num_island << endl;
	cout <<"info.gw_alpha        = " << info.gw_alpha << endl;
	cout <<"info.gw_beta         = " << info.gw_beta << endl;
	cout <<"info.lw_alpha        = " << info.lw_alpha << endl;
	cout <<"info.lw_beta         = " << info.lw_beta << endl;
	cout <<"info.bgh_alpha       = " << info.bgh_alpha << endl;
	cout <<"info.bgh_beta        = " << info.bgh_beta << endl;
	cout <<"info.bgl_alpha       = " << info.bgl_alpha << endl;
	cout <<"info.bgl_beta        = " << info.bgl_beta << endl;
	cout <<"info.prob_cross      = " << info.prob_cross << endl;
	cout <<"info.prob_exchange   = " << info.prob_exchange << endl;
	cout <<"info.prob_mutation   = " << info.prob_mutation << endl;
	cout <<"info.beta_log_upper  = " << info.beta_log_upper << endl;
	cout <<"info.beta_log_lower  = " << info.beta_log_lower << endl;
	// cout <<"info.bg_stddev = " << info.bg_stddev << endl;
	// cout <<"info.shift_stddev = " << info.shift_stddev << endl;
	cout <<"info.intensity_alpha = " << info.intensity_alpha << endl;
	cout <<"info.intensity_beta  = " << info.intensity_beta << endl;
	cout <<"info.intensity_max   = " << info.intensity_max << endl;
	cout <<"info.noiseType       = " << info.noiseType << endl;

    return true;
}

int SetupDataObj3()
{
    d.Vg = newFuncTbl();     // Initialize a table of Voigt function
	// string voigt_ref_file = "../../../../src/VoigtDB.csv";
	// string voigt_ref_file = "./VoigtDB.csv";
	string voigt_ref_file = "VoigtDB.csv";
    LoadFtbl(voigt_ref_file.c_str(), d.Vg);
	return 0;
}

int execDistributedGeneticAlgorithm(int &random_seed, SolutionStruct &solution){
	char fbuf[1024];
    // unsigned int NumOfCores=thread::hardware_concurrency();
	// cout << "the number of core : " << NumOfCores << endl;
	cout << "Genetic Algorithm" << endl;
	// reading config file
	InfoGenetic info;
    if(!readConfig(info)){ cout << "ERROR [config file]" << endl; return -1; }
	info.seed = random_seed;
	// showInfo(info);

	info.p_upper = solution.series.data.x.front();
	info.p_lower = solution.series.data.x.back();
	int num_peaks = solution.series.fitting_parameters.peaks_size;

	// if( ( NumOfCores - 1 ) <  info.num_island){ info.num_island = (int)( NumOfCores / 2 ); }
	// if(info.num_island%2 != 0){ cout << "the number of island is add!" << endl; return 0; }
	if(info.num_island%2 != 0){ cout << "Error! Number of island must be even!" << endl; return 0; }

	// // int NumOfCoresLimit = NumOfCores/2;
	// // int NumOfCoresLimit = 1;
	// int NumOfCoresLimit = NumOfCores - 1;
	// NumOfCoresLimit = max(1, NumOfCoresLimit);
	// printf("Total number of cores = %d\n", NumOfCores);
	// printf("Upper limit of number of cores = %d\n", NumOfCoresLimit);
	// printf("Number of island = %d\n", info.num_island);

	int numThreadLimit = get_numThreadLimit();

    int num_island = info.num_island;
	vector< InputClass> inputs(num_island);
	// DataObj d = SetupDataObj();

    for (int num=0;num<num_island;num++){

		vector<SolutionStruct> geneticParames;
		geneticParames.resize(info.individual);
		for(int i=0;i<info.individual;i++){
			geneticParames[i] = solution;
		}
		inputs[num].geneticParames = geneticParames;
		inputs[num].info = info;
        inputs[num].info.d = d;
		inputs[num].info.beta = pow( 10.0, random_mt(info.seed, info.beta_log_upper, info.beta_log_lower) );
		inputs[num].best_energy = 1e300;
		inputs[num].log.reserve( info.MaxOfLoop * info.generation );
		
		initializeParameters(inputs[num].geneticParames, inputs[num].info);

	}

    // vector<thread> thread_list(num_island);
    // vector<thread> thread_list;
    for(int loop=0;loop<info.MaxOfLoop;loop++){

		cout << "[ " << (loop+1) << " ] -----------------------------------" << endl;

		// Orginal loop
        // for(int i=0;i<num_island;i++){ thread_list[i] = thread( execGeneticAlgorithm, std::ref(inputs[i]) ); }
        // for(thread &th : thread_list){ th.join(); }

		// Serial loop
		// for(int i=0;i<num_island;i++){ execGeneticAlgorithm(inputs[i]); }

		// Partially parallel loop
		int j = 0;
		do {
			vector<thread> thread_list;
			// vector<thread> thread_list(NumOfCoresLimit);
			for(int i=0; i < numThreadLimit; i++){
				// printf("dbg. parallel loop for j = %d, core ID = %d\n", j, i);
				thread_list.emplace_back(thread( execGeneticAlgorithm, std::ref(inputs[j]) ));
				// thread_list[i] = thread( execGeneticAlgorithm, std::ref(inputs[j]) );
				j++;
				if (j >= num_island){break;}
			}
			for(thread &th : thread_list){ th.join(); }
		} while(j < num_island - 1);

		// printf("dbg. finish loop\nPress enter\n"); getchar();

        vector<int> v = makeRandomIndexList(num_island, info.seed);

		for(int i=0;i<num_island;i=i+2) {
			for(int j=0;j<info.individual;j++) {
				if(random_mt(info.seed, 1.0, 0.0)<info.prob_exchange) {
					swap(inputs[v[i]].geneticParames[j], inputs[v[i+1]].geneticParames[j]);
				}
			}
		}

		if( ( loop + 1 ) % 5 == 0 && loop != 0 ){
			std::ostringstream ss;
			ss << "result_"<< std::setw(4) << std::setfill('0') << (loop + 1);
			//std::filesystem::create_directory(ss.str());
			//std::experimental::filesystem::create_directory(ss.str());


			// Make folder "result_####" without using <filesystem>
			sprintf(fbuf, "mkdir -p %s\n", ss.str().c_str());
			printf("%s", fbuf);
			if(system(fbuf)){
				printf("Error: mkdir failed. %s\n", fbuf);
			}

			SolutionStruct tmp_result;
			calcMinimumEnergyResult(num_island, tmp_result, inputs);
			saveFiles(ss.str(), tmp_result, inputs, d);
		}

    }

	saveAllIsland(inputs);

    SolutionStruct result;
	calcMinimumEnergyResult(num_island, result, inputs);

	cout << "< save files >" << endl;
	cout << "loss = " << result.energy << endl;
	string s = "result";
	//std::filesystem::create_directory(s);
	//std::experimental::filesystem::create_directory(s);

	// Make folder "result" without using <filesystem>
	sprintf(fbuf, "mkdir -p %s\n", s.c_str());
	printf("%s", fbuf);
	if(system(fbuf)){
		printf("Error: mkdir failed. %s\n", fbuf);
	}


    saveFiles(s, result, inputs, d);

    return 0;
}

int initializeParameters(vector<SolutionStruct> &geneticParames, InfoGenetic &info){
	for (int i=0;i<info.individual;i++){
		SamplingFromPriorDistributions(geneticParames[i], info);
	}

    return 0;
}

int SamplingFromPriorDistributions(SolutionStruct &geneticParame, InfoGenetic &info){
	double energy_step = abs( geneticParame.series.data.x[0] - geneticParame.series.data.x[1] );

	int num_peaks = geneticParame.series.fitting_parameters.peaks_size;

	geneticParame.series.fitting_parameters.parameters.resize(num_peaks);
    for(int k = 0; k < num_peaks; k++){
		//geneticParame.series.fitting_parameters.parameters[k].area = random_gamma_mt(info.seed, info.intensity_alpha, info.intensity_beta);
		double intensity;
		do {
			intensity = random_gamma_mt(info.seed, info.intensity_alpha, info.intensity_beta);
			//cout << "intensity (random)" << intensity << endl; 
		} while (intensity > info.intensity_max);
		geneticParame.series.fitting_parameters.parameters[k].area = intensity;
		//cout << "intensity (sampled)" << geneticParame.series.fitting_parameters.parameters[k].area << endl; 

		geneticParame.series.fitting_parameters.parameters[k].position = random_mt(info.seed, info.p_upper, info.p_lower);
		geneticParame.series.fitting_parameters.parameters[k].gauss_width = max( random_gamma_mt(info.seed, info.gw_alpha, info.gw_beta), 1e-8 );
		geneticParame.series.fitting_parameters.parameters[k].lorentz_width = max( random_gamma_mt(info.seed, info.lw_alpha, info.lw_beta), 1e-8 );
    }

	
	sortPositionParameter(geneticParame.series.fitting_parameters.parameters);

	// double bgh_mean = geneticParame.series.data.y.front();
	// double bgl_mean = geneticParame.series.data.y.back();
	// geneticParame.series.high = max( random_gauss_mt(info.seed, bgh_mean, info.bg_stddev), 1e-8);
	// geneticParame.series.low = max( random_gauss_mt(info.seed, bgl_mean, info.bg_stddev), 1e-8);

	
	

	// double highMax = (info.bgh_alpha - 1.0) * info.bgh_beta; // mode of gamma distribution
	double highMax = info.bgh_alpha * info.bgh_beta; // mean of gamma distribution
	double high;
	do {
		high = random_gamma_mt(info.seed, info.bgh_alpha, info.bgh_beta);
		// cout << "high (random)" << high << endl; 
	} while (high > highMax);
	geneticParame.series.high = high;

	// double lowMax = (info.bgl_alpha - 1.0) * info.bgl_beta; // mode of gamma distribution
	double lowMax = info.bgl_alpha * info.bgl_beta; // mean of gamma distribution
	double low;
	do {
		low = random_gamma_mt(info.seed, info.bgl_alpha, info.bgl_beta);
		// cout << "low (random)" << low << endl; 
	} while (low > lowMax);
	geneticParame.series.low = low;

	// cout << "dbg" << endl; getchar();

	return 0;
}

int sortPositionParameter(vector<CVP> &parameters){
	int num_peaks = parameters.size();

	vector<double> position(num_peaks);
    for(int k = 0; k < num_peaks; k++){
		position[k] = parameters[k].position;
	}
	std::sort( position.begin(), position.end() );
	for(int k = 0; k < num_peaks; k++){
		parameters[k].position = position[k];
	}

	return 0;
}

int execGeneticAlgorithm(InputClass& input){

	vector<SolutionStruct> geneticParames = input.geneticParames;
	InfoGenetic info = input.info;
	vector<string> labels = input.labels;
	vector<double> energy_list(info.individual);
	size_t min_energy_id;
	stringstream ss;

	for(int g=0;g<info.generation;g++){

		calcFitness(geneticParames, info, energy_list);
		
		double print_loss = ( *min_element(energy_list.begin(), energy_list.end()) );
		ss << g << "," << print_loss << endl;

		min_energy_id = distance(
			energy_list.begin(),
			min_element(energy_list.begin(), energy_list.end())
		);
		if( input.best_energy > energy_list[min_energy_id] ){
			input.bestGeneticParame = geneticParames[min_energy_id];
			input.best_energy = energy_list[min_energy_id];
		}

		geneticParames[0] = input.bestGeneticParame;

		geneticParames = selectIndividual(geneticParames, info);
		geneticParames = crossOverGene(geneticParames, info);
		geneticParames = mutateGene(geneticParames, info);

		// Progress
		// if( (g + 1)%200 == 0 || g == 0){
		// 	cout <<"ID:" << std::this_thread::get_id()
		// 	<< ", loop=" << (g + 1) << ", loss=" << print_loss << endl;
		// }

		input.log.push_back(input.best_energy);

	}

	input.geneticParames = geneticParames;

	return 0;
}

int calcFitness(vector<SolutionStruct> &geneticParames, InfoGenetic info, vector<double> &energy_list)
{
	for(int ind=0;ind<info.individual;ind++){
		geneticParames[ind].calc_energy(info);
		energy_list[ind] = geneticParames[ind].energy;
	}

	return 0;
}

vector<SolutionStruct> selectIndividual(vector<SolutionStruct> &geneticParames, InfoGenetic &info){

	vector<int> index_list;
	vector<vector<int> > local_index_list;
	calcPorbList(geneticParames, info, index_list, local_index_list);

	vector<SolutionStruct> new_geneticParames;
	new_geneticParames.resize(info.individual);

	for(int i = 0; i < info.individual; i++){
		SolutionStruct solution = geneticParames[index_list[i]];
		new_geneticParames[i] = solution;
	}

	return new_geneticParames;

}

vector<double> fitness2prob(vector<double> fitness){
	double sum = accumulate(fitness.begin(), fitness.end(), 0.0);
	for(auto& f : fitness){ f = f / sum; }
	return fitness;
}

vector<int> prob2indexList(const InfoGenetic &info, vector<double> prob){
	int size = prob.size();
	vector<int> index_list(size);
	int index;
	double sum_prob, r;
	for (int i=0;i<size;i++){
		r = random_mt(info.seed, 1.0, 0.0);
		sum_prob = 0.0;
		index = 0;  
		for (auto p : prob){
			sum_prob += p;
			if(r <= sum_prob){ break; }
			else{ index++; }
		}
		index_list[i] = index;
	}

	return index_list;
}

int calcSelectionPressure(InfoGenetic &info, vector<double> &fitness_list){
	double base = 1.0;
	double max_fitness = *max_element(fitness_list.begin(), fitness_list.end());
	double min_fitness = *min_element(fitness_list.begin(), fitness_list.end());
	double scale = 1.0 / (max_fitness - min_fitness);
	for (auto& f : fitness_list){
        // f = scale * ( (info.beta - base) * f + (base * max_fitness - info.beta * min_fitness) );
		f = pow(f, info.beta);
	}

	return 0;
}

int calcPorbList(vector<SolutionStruct> &geneticParames,  InfoGenetic &info,
vector<int> &index_list, vector<vector<int> > &local_index_list){

	vector<double> fitness_list(info.individual), prob_list(info.individual);
	index_list.resize(info.individual);

	for (int i = 0; i < info.individual; i++){
		fitness_list[i] = ( 1.0 / geneticParames[i].energy );
	}

	calcSelectionPressure( info, fitness_list );
	prob_list = fitness2prob( fitness_list );
	index_list = prob2indexList( info, prob_list );

	return 0;
}

vector<SolutionStruct> mutateGene(vector<SolutionStruct> &geneticParames, InfoGenetic &info){
	
	double ran_muta;
	int cnt=0;
	for(int i=0;i<info.individual;i++){
		ran_muta = random_mt(info.seed, 1.0, 0.0);
		if(ran_muta < info.prob_mutation){
			SamplingFromPriorDistributions(geneticParames[i], info);
		}
	}

	return geneticParames;
}

vector<SolutionStruct> crossOverGene(vector<SolutionStruct> &geneticParames, const InfoGenetic &info)
{

	vector<SolutionStruct> new_geneticParames;
	new_geneticParames.resize(info.individual);
	for(int i=0;i<info.individual;i++){
        new_geneticParames[i] = geneticParames[i];
	}

    for(int i = 0; i < info.individual; i += 2 ){
		if( random_mt(info.seed, 1.0, 0.0) < info.prob_cross ){
            vector<SolutionStruct> children = calcCrossover(geneticParames[i], geneticParames[i+1], info);
            for( int k = 0; k < 2; k++ ){ new_geneticParames[i+k] = children[k]; }
        }
    }

    return new_geneticParames;
}

vector<SolutionStruct> calcCrossover(SolutionStruct &struct_dad, SolutionStruct &struct_mam, const InfoGenetic &info){
    vector<SolutionStruct> children(2);
    children[0] = struct_dad;
    children[1] = struct_mam;

	int num_peaks = struct_dad.series.fitting_parameters.peaks_size;
	double energy_step = abs( struct_dad.series.data.x[0] - struct_dad.series.data.x[1] );
	double n = 0.01;

    for( int m=0; m < 2; m++ ){
		vector<CVP> param(num_peaks);
        for(int k = 0; k < num_peaks; k++){
			param[k].area = generateChild(
                struct_dad.series.fitting_parameters.parameters[k].area,
                struct_mam.series.fitting_parameters.parameters[k].area,
                info );
			param[k].area += ( param[k].area * random_mt(info.seed, n, - n) );
			param[k].area = max( param[k].area, 1e-8);

            param[k].position = generateChild(
                struct_dad.series.fitting_parameters.parameters[k].position,
                struct_mam.series.fitting_parameters.parameters[k].position,
                info );
			param[k].position += random_mt(info.seed, 0.05, -0.05);
			// Limit the peak position within [p_lower to p_upper]
			param[k].position = max(min(param[k].position, info.p_upper), info.p_lower);

            param[k].gauss_width = generateChild(
                    struct_dad.series.fitting_parameters.parameters[k].gauss_width,
                    struct_mam.series.fitting_parameters.parameters[k].gauss_width,
                    info );
			param[k].gauss_width += ( param[k].gauss_width * random_mt(info.seed, n, - n) );
            param[k].gauss_width = max( param[k].gauss_width, 1e-8);

			param[k].lorentz_width = generateChild(
                    struct_dad.series.fitting_parameters.parameters[k].lorentz_width,
                    struct_mam.series.fitting_parameters.parameters[k].lorentz_width,
                    info );
			param[k].lorentz_width += ( param[k].lorentz_width * random_mt(info.seed, n, - n) );
            param[k].lorentz_width = max( param[k].lorentz_width, 1e-8);
            
        }
		sortPositionParameter(param);
		children[m].series.fitting_parameters.parameters = param;

		children[m].series.high = generateChild(
			struct_dad.series.high,
			struct_mam.series.high,
			info );
		children[m].series.high = max( children[m].series.high, 1e-8);

		children[m].series.low = generateChild(
			struct_dad.series.low,
			struct_mam.series.low,
			info );
		children[m].series.low = max( children[m].series.low, 1e-8);
    }

    return children;
}

double generateChild(double &dad, double &mam, const InfoGenetic &info){
    double r = random_mt(info.seed, 1.0,0.0);
    double val = r * dad + (1.0 - r) * mam;

    return val;
}

inline double random_mt(const int &random_seed, double upper, double lower) {
	static mt19937 mt(random_seed);
	double r = ((double)mt() + 1.0) / ((double)mt.max() + 2.0);
    return (r * (upper - lower) + lower);
}

inline double random_gauss_mt(const int &random_seed, double mean, double stddev) {
	static mt19937 engine(random_seed);
    normal_distribution<> dist(mean, stddev);
	return dist(engine);
}

inline double random_gamma_mt(const int &random_seed, double alpha, double beta) {
	static mt19937 engine(random_seed);
	gamma_distribution<> dist(alpha, beta);
	return dist(engine);
}

int outputEnergylog(vector<InputClass> &inputs, const string &ofname){
	ofstream ofs;
    stringstream ss;

	ss << "id";
	int loop_size = ( inputs[0].info.MaxOfLoop * inputs[0].info.generation );
	int num_island = inputs[0].info.num_island;
	for(int i=0; i < num_island; i++){
		ss << ',' << "island:" << i;
	}
	ss << endl;
	int cnt = 0;
	for ( int i=0; i < loop_size; i++ ){
		ss << cnt;
		for ( int j=0; j < num_island; j++ ){
			ss << "," << inputs[j].log[i];
		}
		cnt++;
		ss << endl;
	}

	ofs.open(ofname);
	if (ofs) {
		// cout << "Output to " + ofname << endl;
		ofs << ss.str();
		ofs.close();
	}
	else {
		cout << "\n!!! Failed to create a file: " << ofname << " !!!\n" << endl;
	}
	ss.str("");

	return 0;
}

vector<int> makeRandomIndexList(const int &size, const int &seed){
    vector<int> v(size);
    iota(v.begin(), v.end(), 0);
    mt19937 engine(seed);
    shuffle(v.begin(), v.end(), engine);

    return v;
}

int calcMinimumEnergyResult(int &num_island, SolutionStruct &result, vector<InputClass> &inputs){
	// SolutionStruct result;
	double min_energy = 1e300;
	for (int i=0;i<num_island;i++){
		if( min_energy > inputs[i].best_energy ){
			result = inputs[i].bestGeneticParame;
			min_energy = inputs[i].best_energy;
		}
	}

	return 0;
}

int saveFiles(const string ss, SolutionStruct &solution, vector<InputClass> &inputs, DataObj d){

	std::ostringstream ofname;
	ofname << ss << "/" << "opt_result.csv";
	solution.series.outputFitting(ofname.str(), &d);
    solution.series.fitting_parameters.outputCommonParameters(ss + "/" + "opt_parameter.csv");

    outputEnergylog(inputs, ss + "/" + "energy_log.csv");

	return 0;
}

int saveAllIsland(vector<InputClass> &inputs){

	// std::ostringstream ofname;
	// ofname << "./opt_parameters_all_island.csv";
	string ofname = "./opt_parameters_all_island.csv";

    stringstream ss;

    // ss << "compound,shift,gauss_width,lorentz_width,sigma";
    // ss << "compound,height,position,gauss_width,lorentz_width";
    // ss << endl;

	for(auto &in : inputs){
		const Series* s = &in.bestGeneticParame.series;
		ss << s->fitting_parameters.parameters.size() << "," 
		<< s->high << ","
		<< s->low;
		for(auto &p : s->fitting_parameters.parameters){
			ss << "," << p.area 
			<< "," << p.position
			<< "," << p.gauss_width
			<< "," << p.lorentz_width;
		}
		ss << endl;
	}

	// ofname = kwd;
	ofstream ofs;
	ofs.open(ofname);
	if (ofs) {
		// cout << "Output to " + ofname << endl;
		ofs << ss.str();
		ofs.close();
	}
	else {
		cout << "\n!!! Failed to create a file: " << ofname << " !!!\n" << endl;
	}
	ss.str("");
	// cout << "debug20210118" << endl; exit(1);

	return 0;
}

int get_numThreadLimit(){
	// Start optimization using thread
    unsigned int NumOfCores=thread::hardware_concurrency();
	printf("Total number of cores = %d\n", NumOfCores);
	int omp_num_threads = get_omp_num_threads();
	printf("OMP_NUM_THREADS = %d\n", omp_num_threads);

	int numThreadLimit = NumOfCores - 1;
	numThreadLimit = max(1, numThreadLimit);
	if(omp_num_threads > 0){
		numThreadLimit = min(numThreadLimit, omp_num_threads);
	}

	printf("Upper limit of number of cores = %d\n", numThreadLimit);
	return numThreadLimit;
}

int get_omp_num_threads() {
	int num_threads = 0;
	const char *tmp = getenv("OMP_NUM_THREADS");
	string env_var(tmp ? tmp : "");
	if (env_var.empty()) {
		cerr << "[ERROR] No such variable found!" << endl;
		//exit(EXIT_FAILURE);
	} else {
		num_threads = stod(env_var);
	}
	return num_threads;
}

