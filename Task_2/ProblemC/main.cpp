#include <utility>
#include <algorithm>
#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>

using namespace std;


class SuffixArray {
	private:
		vector<unsigned> buildCycle();
		string str;
		vector<char> alphabet;

	public:
		SuffixArray(string str, vector<char> alpha) : str(std::move(str)), alphabet(std::move(alpha)){};
		vector<unsigned> buildSuffixArray ();
		vector<unsigned> buildLCP();
	
		void generateCounts (unordered_map<char, unsigned> &alpha_hash, vector<unsigned> &transpositions, vector<unsigned> &classes) const;
		void makeClasses (vector<unsigned> &transpositions, vector<unsigned> &classes) const;
		vector<unsigned>& mainCycle (vector<unsigned> &transpositions, vector<unsigned> &classes, unsigned classes_cnt) const;
		void reversedModular (const vector<unsigned> &suffix_in_alpha, const vector<unsigned> &suffix_arr, unsigned curr_suff, vector<unsigned> &lcp, unsigned &prefix_len) const;
};

vector<unsigned> SuffixArray::buildCycle() {
	vector<unsigned> transpositions(str.size(), 0);
	vector<unsigned> classes(str.size(), 0);

	unordered_map<char, unsigned> alpha_hash;

	unsigned classes_cnt;
	generateCounts(alpha_hash, transpositions, classes);

	transpositions = mainCycle(transpositions, classes, classes_cnt);

	return transpositions;
}

vector<unsigned> &SuffixArray::mainCycle(vector<unsigned> &transpositions, vector<unsigned> &classes, unsigned classes_cnt) const {
	for (unsigned k = 0; (1 << k) < str.size(); k ++) {
		vector<unsigned> trans(str.size(), 0);

		for (unsigned i = 0; i < str.size(); i ++) {
			if (transpositions[i] < (1 << k))
				trans[i] = str.size() + transpositions[i] - (1 << k);

			else
				trans[i] = transpositions[i] - (1 << k);
		}

		vector<unsigned> count(str.size());

		for (unsigned i = 0; i < str.size(); i ++)
			count[classes[trans[i]]] ++;

		for (unsigned i = 1; i < str.size(); i ++)
			count[i] += count[i - 1];

		for (int i = str.size() - 1; i >= 0; i --)
			transpositions[-- count[classes[trans[i]]]] = trans[i];

		classes_cnt = 1;
		vector<unsigned> new_classes(str.size(), 0);

		for (unsigned i = 1; i < str.size(); i ++) {
			unsigned firstHalf1 = (transpositions[i] + (1 << k)) % str.size();

			unsigned firstHalf2 = (transpositions[i - 1] + (1 << k)) % str.size();

			if (classes[transpositions[i]] != classes[transpositions[i - 1]] || classes[firstHalf1] != classes[firstHalf2])
				classes_cnt ++;

			new_classes[transpositions[i]] = classes_cnt - 1;
		}

		classes = new_classes;
	}
	return transpositions;
}

void SuffixArray::generateCounts(unordered_map<char, unsigned> &alpha_hash, vector<unsigned> &transpositions,
							 vector<unsigned> &classes) const {

	for (unsigned i = 0; i < alphabet.size(); i ++)
		alpha_hash[alphabet[i]] = i;

	vector<unsigned> count(alphabet.size(), 0);

	for (char i : str)
		count[alpha_hash[i]] ++;

	for (unsigned i = 1; i < alphabet.size(); i ++)
		count[i] += count[i - 1];

	for (unsigned i = 0; i < str.size(); i ++)
		transpositions[--count[alpha_hash[str[i]]]] = i;

	makeClasses(transpositions, classes);

}

void SuffixArray::makeClasses(vector<unsigned> &transpositions, vector<unsigned> &classes) const {
	unsigned classesCount = 1;

	for (unsigned i = 1; i < str.size(); i ++) {
		if (str[transpositions[i]] != str[transpositions[i - 1]])
			classesCount ++;

		classes[transpositions[i]] = classesCount - 1;
	}
}


vector<unsigned> SuffixArray::buildSuffixArray() {
	string new_str_with_sharp = str + '#';
	vector<char> new_alpha(alphabet.size() + 1);

	new_alpha[0] = '#';

	for (unsigned i = 0; i < alphabet.size(); i ++)
		new_alpha[i + 1] = alphabet[i];

	vector<unsigned> cyclic = buildCycle();
	vector<unsigned> suffix_arr;

	for (unsigned el: cyclic)
		if (el != new_str_with_sharp.size() - 1)
			suffix_arr.push_back(el);

	return suffix_arr;
}


vector<unsigned> SuffixArray::buildLCP () {
	vector<unsigned> lcp(str.size() - 1);
	vector<unsigned> suffix_in_alpha(str.size());
	vector<unsigned> suffix_arr = buildSuffixArray();

	for (unsigned i = 0; i < str.size(); i ++)
		suffix_in_alpha[suffix_arr[i]] = i;

	unsigned prefix_len = 0;

	for (unsigned curr_suff = 0; curr_suff < str.size(); curr_suff ++) {
		if (prefix_len > 0)
			prefix_len --;

		if (suffix_in_alpha[curr_suff] == str.size() - 1)
			prefix_len = 0;

		else {
			reversedModular(suffix_in_alpha, suffix_arr, curr_suff, lcp, prefix_len);
		}

	}

	return lcp;
}


void SuffixArray::reversedModular(const vector<unsigned> &suffix_in_alpha, const vector<unsigned> &suffix_arr,
							  unsigned curr_suff, vector<unsigned> &lcp, unsigned &prefix_len) const {

	unsigned next_suff = suffix_arr[suffix_in_alpha[curr_suff] + 1];

	while (max<unsigned>(curr_suff, next_suff) + prefix_len < str.size()
				   && str[curr_suff + prefix_len] == str[next_suff + prefix_len])
				prefix_len ++;

	lcp[suffix_in_alpha[curr_suff]] = prefix_len;
}


//--------------------------------------------------------------------------------------------


inline bool doesPrefixOfGivenSuffixThatIncludesInFirstSubstringExists (const unsigned &suffixIndex, const unsigned &strLen) {
	return suffixIndex < strLen + 1;
}

string getKthSubstring (const string &str1, const string &str2, const size_t &k) {
	// Precalc?:
	string goodString = str1 + '#' + str2 + '$';
	const vector<char> alphabet = {'#', '$', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z'};
	SuffixArray suffArr(goodString, alphabet);

	// Building suffix Array:
	vector<unsigned> suffixArray = suffArr.buildSuffixArray();
	vector<unsigned> lcp = suffArr.buildLCP();


	// Init:
	size_t currentCount = 0;
	size_t lcpValue = 0;
	int from = -1;
	int length = -1;


	// Main cycle:
	for (size_t i = 2; i < goodString.size() - 1; i ++) {
		const unsigned index1 = suffixArray[i];
		const unsigned index2 = suffixArray[i + 1];

		// If there is a common prefix of two suffixes, THAT ARE IN THE FIRST STRING, we need to refresh value!:
		if ((str1.size() + 1 > index1) != (str1.size() + 1 > index2)) {	
			// If cur val is greater than the prev one, we need to refresh it to!
			// (And count some common substrings, that were met in
			// interval from lcpValue to curLcp (lcp[i])):
			if (lcp[i] > lcpValue) {
				currentCount += lcp[i] - lcpValue;
			}

			lcpValue = lcp[i];
		}

		else {
			lcpValue = min<unsigned>(lcp[i], lcpValue);
		}

		// If we proskochili k common substrings, we need to break our cycle:
		if (k <= currentCount) {
			from = suffixArray[i];
			length = lcp[i] + k - currentCount;
			break;
		}
	}


	// Returning string:
	if (from == -1)
		return "-1";

	else
		return goodString.substr(from, length);
}

void readData (istream &is, string &str1, string &str2, size_t &k) {
	is >> str1 >> str2 >> k;
}

void writeData (ostream &os, const string &str) {
	os << str << "\n";
}


//--------------------------------------------------------------------------------------------

int main () {
	string str1, str2;
	size_t k;
	readData(cin, str1, str2, k);
	string str = getKthSubstring(str1, str2, k);
	writeData(cout, str);
	return 0;
}
