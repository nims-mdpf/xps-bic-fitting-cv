#include "main.h"

int main(int argc, char *argv[]){
    int random_seed = stoi(argv[1]);
	cout << "random seed = " << random_seed << endl;
	int peak_size = stoi(argv[2]);
	cout << "peak size = " << peak_size << endl;

    Series series;
	cout << "Reading [file] >> " << argv[3] << endl;
    readcsv(argv[3], series.data);

	// SetupDataObj3();

    SolutionStruct solution;
    solution.series = series;
	solution.series.fitting_parameters.peaks_size = peak_size;

    execDistributedGeneticAlgorithm(random_seed, solution);

    return 0;
}

int readcsv(const string &fname, Spectrum &data)
{
	ifstream ifs;
	ifs.open(fname.c_str());
	if (ifs.fail()) {
		return 1;
	}
	data.x.clear();
	data.y.clear();
	string line;

	while (getline(ifs, line)) {
		replace(line.begin(), line.end(), ',', ' ');
		if (line.empty() || line[0] == '#') continue;
		stringstream buf(line);
		double xt = 0.0, yt = 0.0;
		if (buf >> xt >> yt) {
				data.x.push_back(xt);
				data.y.push_back(yt);
		}
		else {
			cout << "Error to read a line= " << buf.str() << endl;
		}
	}
	data.size = data.x.size();

	return 0;
}

int readNodeFile(const string &fname, vector<Node> &node){
    ifstream ifs;
	ifs.open(fname.c_str());
	if (ifs.fail()) { return 1; }
    string line;

    vector<string> label;
	getline(ifs, line);
	while (getline(ifs, line)) {
        line.erase(remove(line.begin(), line.end(),' '), line.end());
		if (line[0] == '#') continue;
		replace(line.begin(), line.end(), ',', ' ');
		stringstream buf(line);

		int a, b;
		buf >> a >> b;
		Node n(a, b);
		node.push_back( n );
	}

    return 0;
}