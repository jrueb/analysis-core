#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <iostream>
#include <iomanip>
#include <fstream>
#include <iterator>
#include <vector>
#include <chrono>
#include <algorithm>

#include "Analysis/Core/interface/Analysis.h"
#include "TFile.h" 
#include "TFileCollection.h"
#include "TChain.h"
#include "TH1.h"

using namespace analysis;
using namespace analysis::tools;

//Maximum number of jets allowed for configuration.
const unsigned int MAX_JETS = 4;
//Print configuration to cout on start.
const bool DUMP_PARAM = true;

//A helper function to simplify printing vectors to cout.
template<class T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v) {
	copy(v.begin(), v.end(), std::ostream_iterator < T > (os, " "));
	return os;
}

//Convert an int to a string in the format nth.
std::string int_to_count(int n) {
	std::string res = std::to_string(n);
	if (res.back() == '1')
		res += "st";
	else if (res.back() == '2')
		res += "nd";
	else if (res.back() == '3')
		res += "rd";
	else
		res += "th";
	return res;
}

//Helper class to record and print statistics on the cuts.
class CutFlow {
public:
	CutFlow(const std::vector<std::string>& cut_descriptions,
			const std::string& title = "Cut Flow") {
		_title = title;
		_cut_descriptions = cut_descriptions;
		_nsel.assign(cut_descriptions.size(), 0);
		_cur_cut = 0;
		started = false;
	}

	//Inform that the next even is being processed.
	void next_event() {
		_cur_cut = 0;
		if (not started) {
			start_time = std::chrono::steady_clock::now();
			started = true;
		}
	}

	//Inform that a cut is done if shouldCut else that the even is
	//being processed further. For convenience returns shouldCut.
	bool do_cut(bool shouldCut) {
		if (not shouldCut)
			++_nsel[_cur_cut];
		++_cur_cut;
		return shouldCut;
	}

	//Inform that all cuts have been done and even was cut after n cuts.
	//Returns true if the event wasn't cut.
	bool do_cut_all(unsigned int n) {
		for (unsigned int i = 0; i < n; ++i)
			++_nsel[i];
		return n != _nsel.size();
	}

	const std::string& get_title() const {
		return _title;
	}

	const std::vector<std::string>& get_descriptions() const {
		return _cut_descriptions;
	}

	const std::vector<unsigned int>& get_nsel() const {
		return _nsel;
	}

	float get_effeciency() {
		if (_nsel.front() == 0)
			return 1.;
		return (float) _nsel.back() / _nsel.front();
	}

	//Return the time in ms since the first event was processed.
	int get_duration() {
		if (not started)
			return 0;
		return std::chrono::duration_cast < std::chrono::milliseconds
				> (std::chrono::steady_clock::now() - start_time).count();
	}

private:
	std::string _title;
	unsigned int _cur_cut;
	std::vector<std::string> _cut_descriptions;
	std::vector<unsigned int> _nsel;
	bool started;
	std::chrono::time_point<std::chrono::steady_clock> start_time;
};
std::ostream& operator<<(std::ostream& os, const CutFlow& obj) {
	size_t sp[] = { 25, 12, 12, 12 };
	const std::vector<std::string>& desc = obj.get_descriptions();
	const std::vector<unsigned int>& nsel = obj.get_nsel();
	sp[0] = std::max(sp[0], obj.get_title().size() + 3);
	for (const auto& d : desc)
		sp[0] = std::max(sp[0], d.size() + 3);

	os << std::setw(sp[0]) << std::left << obj.get_title() << std::setw(sp[1])
			<< std::right << "# events" << std::setw(sp[2]) << "absolute"
			<< std::setw(sp[3]) << "relative";
	for (unsigned int i = 0; i < nsel.size(); ++i) {
		float fracAbs = (float) nsel[i] / nsel[0];
		float fracRel = 1.;
		if (i > 0)
			fracRel = (float) nsel[i] / nsel[i - 1];

		os << std::endl;
		os << std::setw(sp[0]) << std::left << desc[i] << std::setw(sp[1])
				<< std::right << nsel[i] << std::fixed << std::setw(sp[2])
				<< std::setprecision(6) << fracAbs << std::setw(sp[3])
				<< fracRel;

	}
	return os;
}

std::map<std::string, TH1F*> create_histograms(unsigned int njets) {
	std::map<std::string, TH1F*> h1;
	h1["total"] = new TH1F("total", "", 1, 0., 1.);
	h1["totalSelected"] = new TH1F("totalSelected", "", 1, 0., 1.);
	h1["n"] = new TH1F("n", "", 30, 0, 30);
	h1["n_csv"] = new TH1F("n_csv", "", 30, 0, 30);
	h1["n_ptmin20"] = new TH1F("n_ptmin20", "", 30, 0, 30);
	h1["n_ptmin20_csv"] = new TH1F("n_ptmin20_csv", "", 30, 0, 30);
	h1["n_ptmin30_csv"] = new TH1F("n_ptmin30_csv", "", 30, 0, 30);
	for (unsigned int i = 0; i < njets; ++i) {
		if (i < 2) {
			h1[Form("pt_%i", i)] = new TH1F(Form("pt_%i", i), "", 100, 0, 1000);
			h1[Form("eta_%i", i)] = new TH1F(Form("eta_%i", i), "", 100, -5, 5);
			h1[Form("phi_%i", i)] = new TH1F(Form("phi_%i", i), "", 100, -4, 4);
			h1[Form("btag_%i", i)] = new TH1F(Form("btag_%i", i), "", 100, 0,
					1);

			h1[Form("pt_%i_csv", i)] = new TH1F(Form("pt_%i_csv", i), "", 100,
					0, 1000);
			h1[Form("eta_%i_csv", i)] = new TH1F(Form("eta_%i_csv", i), "", 100,
					-5, 5);
			h1[Form("phi_%i_csv", i)] = new TH1F(Form("phi_%i_csv", i), "", 100,
					-4, 4);
			h1[Form("btag_%i_csv", i)] = new TH1F(Form("btag_%i_csv", i), "",
					100, 0, 1);
		} else {
			h1[Form("pt_%i", i)] = new TH1F(Form("pt_%i", i), "", 50, 0, 200);
			h1[Form("eta_%i", i)] = new TH1F(Form("eta_%i", i), "", 100, -5, 5);
			h1[Form("phi_%i", i)] = new TH1F(Form("phi_%i", i), "", 100, -4, 4);
			h1[Form("btag_%i", i)] = new TH1F(Form("btag_%i", i), "", 100, 0,
					1);

			h1[Form("pt_%i_csv", i)] = new TH1F(Form("pt_%i_csv", i), "", 50, 0,
					200);
			h1[Form("eta_%i_csv", i)] = new TH1F(Form("eta_%i_csv", i), "", 100,
					-5, 5);
			h1[Form("phi_%i_csv", i)] = new TH1F(Form("phi_%i_csv", i), "", 100,
					-4, 4);
			h1[Form("btag_%i_csv", i)] = new TH1F(Form("btag_%i_csv", i), "",
					100, 0, 1);
		}
	}
	h1["m12"] = new TH1F("m12", "", 50, 0, 1000);
	h1["m12_csv"] = new TH1F("m12_csv", "", 50, 0, 1000);

	return h1;
}

inline std::vector<Jet*> get_jets_loose(Analysis* analysis) {
	// Jets - std::shared_ptr< Collection<Jet> >
	auto slimmedJets = analysis->collection < Jet > ("Jets");
	std::vector<Jet*> selectedJets;
	for (int j = 0; j < slimmedJets->size(); ++j) {
		Jet* jet = &slimmedJets->at(j);
		if (jet->idLoose())
			selectedJets.push_back(jet);
	}
	return selectedJets;
}

inline bool select_kinematic(const std::vector<Jet*>& jets, unsigned int njets,
		const float* ptmin, const float* etamax) {
	for (unsigned int j = 0; j < njets; ++j) {
		Jet* jet = jets[j];
		if (jet->pt() < ptmin[j] || fabs(jet->eta()) > etamax[j])
			return false;
	}
	return true;
}

inline bool select_deltaR(const std::vector<Jet*>& jets, unsigned int njets,
		float dRmin) {
	for (unsigned int j1 = 0; j1 < njets - 1; ++j1) {
		const Jet& jet1 = *jets[j1];
		for (unsigned int j2 = j1 + 1; j2 < njets; ++j2) {
			const Jet& jet2 = *jets[j2];
			if (jet1.deltaR(jet2) < dRmin)
				return false;
		}
	}
	return true;
}

inline bool select_btag(const std::vector<Jet*>& jets, unsigned int njets,
		bool deepb, bool isbbbb, const float* btagmin, float nonbtag) {
	for (unsigned int j = 0; j < njets; ++j) {
		Jet* jet = jets[j];
		float btag;
		if (deepb)
			btag = jet->btag("btag_deepb") + jet->btag("btag_deepbb");
		else
			btag = jet->btag();
		if ((j < njets - 1 and btag < btagmin[j])
				or (j == njets - 1 and isbbbb and btag < btagmin[j])
				or (j == njets - 1 and not isbbbb and btag > nonbtag)) {
			return false;
		}
	}
	return true;
}

inline unsigned int get_num_matched(const std::vector<Jet*>& jets,
		unsigned int njets, const std::vector<std::string>& triggerObjects) {
	unsigned int n;
	for (unsigned int j = 0; j < njets; ++j) {
		Jet* jet = jets[j];
		n = 0;
		for (const auto& triggerObject : triggerObjects) {
			if (not jet->matched(triggerObject))
				return n;
			++n;
		}
	}
	return n;
}

int main(int argc, char* argv[]) {
	unsigned int njets;
	std::string config_file, input_list, output_file, json_file;
	std::string jetTreePath, triggerResultsPath, triggerBranch;
	std::vector < std::string > triggerObjects;
	bool isbbbb, isMC;
	bool deepb;
	float ptmin[MAX_JETS];
	float btagmin[MAX_JETS];
	float nonbtag;
	float etamax[MAX_JETS];
	float dRmin;
	float detamax;

	try {
		// Declare a group of options that will be 
		// allowed only on command line
		po::options_description generic("Generic options");
		generic.add_options()("help,h", "Produce help message.")("config,c",
				po::value < std::string > (&config_file),
				"Name of a file of a configuration.");

		// Declare a group of options that will be 
		// allowed both on command line and in
		// config file
		po::options_description config("Configuration");
		config.add_options()("output",
				po::value < std::string
						> (&output_file)->default_value("histograms.root"),
				"Name of the ouput root file for histograms.")("json",
				po::value < std::string > (&json_file),
				"Name of the JSON file. If supplied, it is assumed "
						"that input is data and not MC.")("jettreepath",
				po::value < std::string
						> (&jetTreePath)->default_value("MssmHbb/Events/"
								"slimmedJetsPuppiReapplyJEC"),
				"Path within the input root files to the jet tree.")(
				"trigrespath", po::value < std::string > (&triggerResultsPath),
				"Path within the input root files to the trigger "
						"results.")("trigbranch",
				po::value < std::string > (&triggerBranch),
				"Branch within tregrespath to use as a trigger.")("trigobj",
				po::value < std::vector
						< std::string >> (&triggerObjects)->composing(),
				"Trigger object to match against.")("njets",
				po::value<unsigned int>(&njets)->default_value(MAX_JETS),
				"Minimum number of jets")("deepb",
				"Use DeepFlavour btag discriminator value instead"
						" of the normal btag value.");
		for (unsigned int i = 0; i < MAX_JETS; ++i) {
			std::string help_text = "Minimum pt of the " + int_to_count(i + 1)
					+ " leading jet.";
			config.add_options()(
					(std::string("ptmin") + std::to_string(i + 1)).c_str(),
					po::value<float>(&ptmin[i])->default_value(0.),
					help_text.c_str());
		}
		for (unsigned int i = 0; i < MAX_JETS; ++i) {
			std::string help_text = "Minimum btag value of the "
					+ int_to_count(i + 1) + " leading jet.";
			config.add_options()(
					(std::string("btagmin") + std::to_string(i + 1)).c_str(),
					po::value<float>(&btagmin[i])->default_value(0.),
					help_text.c_str());
		}
		config.add_options()("nonbtag", po::value<float>(&nonbtag),
				"Maximum btag value for the njets-th leading jet (e.g. for "
						"CR). If supplied, the last btagmin value is ignored.");
		for (unsigned int i = 0; i < MAX_JETS; ++i) {
			std::string help_text = "Maximum eta of the " + int_to_count(i + 1)
					+ " leading jet.";
			config.add_options()(
					(std::string("etamax") + std::to_string(i + 1)).c_str(),
					po::value<float>(&etamax[i])->default_value(1.),
					help_text.c_str());
		}
		config.add_options()("drmin",
				po::value<float>(&dRmin)->default_value(1.),
				"Minimum delta R between the first njets jets.")("detamax",
				po::value<float>(&detamax)->default_value(6.3),
				"Maximum delta eta between the two leading jets");

		// Hidden options, will be allowed both on command line and
		// in config file, but will not be shown to the user.
		po::options_description hidden("Hidden options");
		hidden.add_options()("input-list",
				po::value < std::string
						> (&input_list)->default_value("rootFileList.txt"),
				"Name of the input file list with path to ntuples.");

		// Add options
		po::options_description cmdline_options;
		cmdline_options.add(generic).add(config).add(hidden);

		po::options_description config_file_options;
		config_file_options.add(config).add(hidden);

		po::options_description visible("Allowed options");
		visible.add(generic).add(config);

		//Make it possible to pass input-list name as
		//positional argument.
		po::positional_options_description p;
		p.add("input-list", -1);

		//Parse arguments from command line
		po::variables_map vm;
		po::store(
				po::command_line_parser(argc, argv).options(cmdline_options).positional(
						p).run(), vm);

		if (vm.count("help")) {
			//Print help text.
			std::cout << visible << std::endl;
			return 0;
		}
		//Throw exceptions if applicable
		po::notify(vm);

		if (vm.count("config")) {
			//Parse config file
			std::ifstream ifs(config_file.c_str());
			if (not ifs) {
				std::cout << "Can not open config file: " << config_file
						<< std::endl;
				return 1;
			} else {
				store(parse_config_file(ifs, config_file_options), vm);
				notify(vm);
			}
		}

		//Set booleans if particular arguments are supplied
		isMC = vm.count("json") == 0;
		isbbbb = vm.count("nonbtag") == 0;
		deepb = vm.count("deepb") != 0;

		//Validate configuration
		if (not (njets <= MAX_JETS)) {
			std::cerr << "Invalid value for njets: Outside of interval"
					" [0;" << MAX_JETS << "]." << std::endl;
			return 1;
		}
		for (unsigned int i = 0; i < MAX_JETS; ++i) {
			if (not (0. <= ptmin[i])) {
				std::cerr << "Invalid value for ptmin" << i + 1
						<< ": Less than 0." << std::endl;
				return 1;
			}
			if (not (0. <= btagmin[i] and btagmin[i] <= 1.)) {
				std::cerr << "Invalid value for btagmin" << i + 1
						<< ": Outside of interval [0;1]." << std::endl;
				return 1;
			}
			//TODO: Limits for eta?

		}
		if (vm.count("nonbtag")) {
			if (not (0. <= nonbtag and nonbtag <= 1.)) {
				std::cerr << "Invalid value for nonbtag: "
						"Outside of interval [0;1]." << std::endl;
				return 1;
			}
		}
		if (not (0. <= dRmin)) {
			std::cerr << "Invalid value for drmin: Less than 0." << std::endl;
			return 1;
		}
		if (not (0. <= detamax)) {
			std::cerr << "Invalid value for deltamax: Less than 0."
					<< std::endl;
			return 1;
		}
		std::ifstream input_list_s(input_list.c_str());
		if (not input_list_s) {
			std::cerr << "Can not open input list file: " << input_list
					<< std::endl;
			return 1;
		}
		input_list_s.close();
		if (vm.count("json")) {
			std::ifstream json_file_s(json_file.c_str());
			if (not json_file_s) {
				std::cerr << "Can not open json file: " << json_file
						<< std::endl;
				return 1;
			}
			json_file_s.close();
		}
		std::ofstream output_file_s(output_file.c_str());
		if (not output_file_s) {
			std::cerr << "Can not open output list file: " << output_file
					<< std::endl;
			return 1;
		}
		output_file_s.close();

		if (DUMP_PARAM) {
			//Print all parameters to cout
			std::cout << "njets=" << njets << std::endl;
			if (deepb)
				std::cout << "Using btag DeepFlavour." << std::endl;
			else
				std::cout << "Using normal b-tagging." << std::endl;
			std::cout << "ptmin=" << std::vector<float>(ptmin, ptmin + njets)
					<< std::endl;
			std::cout << "btagmin="
					<< std::vector<float>(btagmin, btagmin + njets)
					<< std::endl;
			if (vm.count("nonbtag"))
				std::cout << "nonbtag=" << nonbtag << std::endl;
			else
				std::cout << "No nonbtag." << std::endl;
			std::cout << "etamax=" << std::vector<float>(etamax, etamax + njets)
					<< std::endl;
			std::cout << "dRmin=" << dRmin << std::endl;
			std::cout << "detamax=" << detamax << std::endl;
			std::cout << "Input list is " << input_list << std::endl;
			std::cout << "Output file is " << output_file << std::endl;
			if (vm.count("json"))
				std::cout << "JSON file is " << json_file << std::endl;
			else
				std::cout << "Running without JSON file." << std::endl;
			std::cout << "jetTreePath=" << jetTreePath << std::endl;
			std::cout << "triggerResultsPath=" << triggerResultsPath
					<< std::endl;
			std::cout << "triggerBranch=" << triggerBranch << std::endl;
			std::cout << "triggerObjects (n=" << triggerObjects.size() << "): "
					<< triggerObjects << std::endl;

		}
	} catch (std::exception& e) {
		std::cerr << e.what() << std::endl;
		return 1;
	}

	TH1::SetDefaultSumw2(); // proper treatment of errors when scaling histograms
	Analysis analysis(input_list);
	analysis.addTree < Jet > ("Jets", jetTreePath);
	analysis.triggerResults(triggerResultsPath);
	for (const auto& obj : triggerObjects)
		analysis.addTree < TriggerObject > (obj, obj.c_str());
	if (not isMC)
		analysis.processJsonFile(json_file);

	std::map<std::string, TH1F*> h1 = create_histograms(njets);

	//Make a list of trigger object names without preceding path
	std::vector < std::string > trigNames;
	for (const auto& obj : triggerObjects) {
		size_t p = obj.rfind('/');
		trigNames.push_back(obj.substr(p == std::string::npos ? 0 : p + 1));
	}

	CutFlow cf(
			{ "JSON", "Triggered", "Triple idloose-jet",
					"Triple jet kinematics", "Delta R(i;j)", "Delta eta(j1;j2)",
					"btagged (" + std::string(njets - 1, 'b')
							+ (isbbbb ? "b)" : "nb)"), "Matched to online j1;j2" });
	CutFlow cf_trig(trigNames, "Trigger Objects");

	std::cout << "This analysis has " << analysis.size() << " events."
			<< std::endl;

	for (int i = 0; i < analysis.size(); ++i) {
		if (i > 0 && i % 100000 == 0)
			std::cout << i << " events processed!" << std::endl;
		analysis.event(i);
		cf.next_event();

		h1["total"]->Fill(0.5);

		//Select only good JSON
		if (cf.do_cut(not isMC and not analysis.selectJson()))
			continue;

		//Trigger selection
		if (cf.do_cut(not analysis.triggerResult(triggerBranch)))
			continue;

		//Require minimum of njets loose jets
		std::vector<Jet*> selectedJets = get_jets_loose(&analysis);
		if (cf.do_cut(selectedJets.size() < njets))
			continue;

		//Fill histrograms before further cuts
		h1["n"]->Fill(selectedJets.size());
		h1["n_ptmin20"]->Fill(
				std::count_if(selectedJets.begin(), selectedJets.end(),
						[](Jet* jet) {return jet->pt() >= 20.;}));
		h1["m12"]->Fill((selectedJets[0]->p4() + selectedJets[1]->p4()).M());
		for (unsigned int j = 0; j < njets; ++j) {
			Jet* jet = selectedJets[j];
			h1[Form("pt_%i", j)]->Fill(jet->pt());
			h1[Form("eta_%i", j)]->Fill(jet->eta());
			h1[Form("phi_%i", j)]->Fill(jet->phi());
			h1[Form("btag_%i", j)]->Fill(jet->btag());
		}

		// Kinematic selection
		if (cf.do_cut(not select_kinematic(selectedJets, njets, ptmin, etamax)))
			continue;

		// Delta R selection
		if (cf.do_cut(not select_deltaR(selectedJets, njets, dRmin)))
			continue;

		// Delta eta selection - 2 leading jets
		if (cf.do_cut(
				fabs(selectedJets[0]->eta() - selectedJets[1]->eta())
						> detamax))
			continue;

		// Btag selction
		if (cf.do_cut(
				not select_btag(selectedJets, njets, deepb, isbbbb, btagmin,
						nonbtag)))
			continue;

		// Match offline to online
		analysis.match<Jet, TriggerObject>("Jets", triggerObjects, 0.5);
		// Are the TWO leading jets matched?
		if (cf.do_cut(
				cf_trig.do_cut_all(
						get_num_matched(selectedJets, 2, triggerObjects))))
			continue;

		// Fill histograms of passed btagging selection
		h1["totalSelected"]->Fill(0.5);
		h1["n_csv"]->Fill(selectedJets.size());
		h1["n_ptmin20_csv"]->Fill(
				std::count_if(selectedJets.begin(), selectedJets.end(),
						[](Jet* jet) {return jet->pt() >= 20.;}));
		h1["n_ptmin30_csv"]->Fill(
				std::count_if(selectedJets.begin(), selectedJets.end(),
						[](Jet* jet) {return jet->pt() >= 30.;}));
		for (unsigned int j = 0; j < njets; ++j) {
			Jet* jet = selectedJets[j];
			h1[Form("pt_%i_csv", j)]->Fill(jet->pt());
			h1[Form("eta_%i_csv", j)]->Fill(jet->eta());
			h1[Form("phi_%i_csv", j)]->Fill(jet->phi());
			h1[Form("btag_%i_csv", j)]->Fill(jet->btag());
		}
		h1["m12_csv"]->Fill(
				(selectedJets[0]->p4() + selectedJets[1]->p4()).M());
	}

	//Print statistics
	std::cout << std::endl;
	std::cout << "Analysis took " << cf.get_duration() << " ms." << std::endl;

	TFile hout(output_file.c_str(), "recreate");
	for (const auto& ih1 : h1) {
		ih1.second->Write();
	}
	hout.Close();

	std::cout << cf << std::endl << std::endl << cf_trig << std::endl;

	// Efficiency
	//TODO: Efficiency and JSON?
	std::cout << std::endl;
	std::cout << "Efficiency: " << cf.get_effeciency() << std::endl;

	return 0;
}
