#include <utility>
#include <algorithm>
#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>

using namespace std;


//--------------------------------------------------------------------------------------------


class SuffArr {
private:
    vector<unsigned> buildCycle();
    string m_sString;
    vector<char> m_vcAlphabet;
public:
    SuffArr(string str, vector<char> alpha) : m_sString(std::move(str)), m_vcAlphabet(std::move(alpha)){};
    vector<unsigned> buildSuffixArray ();
    vector<unsigned> buildLCP();

    void generateCounts(unordered_map<char, unsigned int> &alpha_hash, vector<unsigned int> &transpositions,
                        vector<unsigned int> &classes) const;

    void makeClasses(vector<unsigned int> &transpositions, vector<unsigned int> &classes) const;

    vector<unsigned int> &
    mainCycle(vector<unsigned int> &transpositions, vector<unsigned int> &classes, unsigned int classes_cnt) const;

    void reversedModular(const vector<unsigned int> &suffix_in_alpha, const vector<unsigned int> &suffix_arr,
                         unsigned int curr_suff, vector<unsigned int> &lcp, unsigned int &prefix_len) const;

    void generateNewCycles(unsigned int &classes_cnt, unsigned int &k, vector<unsigned int> &trans,
                           vector<unsigned int> &count, vector<unsigned int> &transpositions,
                           vector<unsigned int> &classes, vector<unsigned int> &new_classes) const;
};


vector<unsigned> SuffArr::buildCycle() {
    vector<unsigned> transpositions(m_sString.size(), 0);
    vector<unsigned> classes(m_sString.size(), 0);

    unordered_map<char, unsigned> alpha_hash;

    unsigned int classes_cnt;
    generateCounts(alpha_hash, transpositions, classes);

    transpositions = mainCycle(transpositions, classes, classes_cnt);

    return transpositions;
}


vector<unsigned int> &SuffArr::mainCycle(vector<unsigned int> &transpositions, vector<unsigned int> &classes,
                                         unsigned int classes_cnt) const {
    for (unsigned k = 0; (1 << k) < m_sString.size(); k ++) {
        vector<unsigned> trans(m_sString.size(), 0);

        for (unsigned i = 0; i < m_sString.size(); i ++) {
            if (transpositions[i] < (1 << k))
                trans[i] = m_sString.size() + transpositions[i] - (1 << k);

            else
                trans[i] = transpositions[i] - (1 << k);
        }

        vector<unsigned> count(m_sString.size());

        vector<unsigned int> new_classes;

        generateNewCycles(classes_cnt, k, trans, count, transpositions, classes, new_classes);

        classes = new_classes;
    }
    return transpositions;
}


void SuffArr::generateNewCycles(unsigned int &classes_cnt, unsigned int &k, vector<unsigned int> &trans,
                                vector<unsigned int> &count, vector<unsigned int> &transpositions,
                                vector<unsigned int> &classes, vector<unsigned int> &new_classes) const {
    for (unsigned i = 0; i < m_sString.size(); i ++)
            count[classes[trans[i]]] ++;

    for (unsigned i = 1; i < m_sString.size(); i ++)
            count[i] += count[i - 1];

    for (int i = m_sString.size() - 1; i >= 0; i --)
            transpositions[-- count[classes[trans[i]]]] = trans[i];

    classes_cnt = 1;
    for (unsigned i = 1; i < m_sString.size(); i ++) {
            unsigned firstHalf1 = (transpositions[i] + (1 << k)) % m_sString.size();

            unsigned firstHalf2 = (transpositions[i - 1] + (1 << k)) % m_sString.size();

            if (classes[transpositions[i]] != classes[transpositions[i - 1]] || classes[firstHalf1] != classes[firstHalf2])
                classes_cnt ++;

            new_classes[transpositions[i]] = classes_cnt - 1;
        }
}


void SuffArr::generateCounts(unordered_map<char, unsigned int> &alpha_hash, vector<unsigned int> &transpositions,
                             vector<unsigned int> &classes) const {

    for (unsigned i = 0; i < m_vcAlphabet.size(); i ++)
        alpha_hash[m_vcAlphabet[i]] = i;

    vector<unsigned> count(m_vcAlphabet.size(), 0);

    for (char i : m_sString)
        count[alpha_hash[i]] ++;

    for (unsigned i = 1; i < m_vcAlphabet.size(); i ++)
        count[i] += count[i - 1];

    for (unsigned i = 0; i < m_sString.size(); i ++)
        transpositions[--count[alpha_hash[m_sString[i]]]] = i;

    makeClasses(transpositions, classes);

}


void SuffArr::makeClasses(vector<unsigned int> &transpositions, vector<unsigned int> &classes) const {
    unsigned classesCount = 1;

    for (unsigned i = 1; i < m_sString.size(); i ++) {
        if (m_sString[transpositions[i]] != m_sString[transpositions[i - 1]])
            classesCount ++;

        classes[transpositions[i]] = classesCount - 1;
    }
}


vector<unsigned> SuffArr::buildSuffixArray() {
    string new_str_with_sharp = m_sString + '#';
    vector<char> new_alpha(m_vcAlphabet.size() + 1);

    new_alpha[0] = '#';

    for (unsigned i = 0; i < m_vcAlphabet.size(); i ++)
        new_alpha[i + 1] = m_vcAlphabet[i];

    vector<unsigned> cyclic = buildCycle();
    vector<unsigned> suffix_arr;

    for (unsigned el: cyclic)
        if (el != new_str_with_sharp.size() - 1)
            suffix_arr.push_back(el);

    return suffix_arr;
}


vector<unsigned> SuffArr::buildLCP () {
    vector<unsigned> lcp(m_sString.size() - 1);
    vector<unsigned> suffix_in_alpha(m_sString.size());
    vector<unsigned> suffix_arr = buildSuffixArray();

    for (unsigned i = 0; i < m_sString.size(); i ++)
        suffix_in_alpha[suffix_arr[i]] = i;

    unsigned prefix_len = 0;

    for (unsigned curr_suff = 0; curr_suff < m_sString.size(); curr_suff ++) {
        if (prefix_len > 0)
            prefix_len --;

        if (suffix_in_alpha[curr_suff] == m_sString.size() - 1)
            prefix_len = 0;

        else {
            reversedModular(suffix_in_alpha, suffix_arr, curr_suff, lcp, prefix_len);
        }

    }

    return lcp;
}


void SuffArr::reversedModular(const vector<unsigned int> &suffix_in_alpha, const vector<unsigned int> &suffix_arr,
                              unsigned int curr_suff, vector<unsigned int> &lcp, unsigned int &prefix_len) const {

    unsigned next_suff = suffix_arr[suffix_in_alpha[curr_suff] + 1];

    while (max<unsigned>(curr_suff, next_suff) + prefix_len < m_sString.size()
                   && m_sString[curr_suff + prefix_len] == m_sString[next_suff + prefix_len])
                prefix_len ++;

    lcp[suffix_in_alpha[curr_suff]] = prefix_len;
}


//--------------------------------------------------------------------------------------------


class Solve {
public:
    Solve();
    void printAns();
private:
    unsigned countDifferentSubstrings(string &str, vector<char> &alphabet);
    void readStr();
    string m_sString;
    vector<char> m_vcAlphabet =
                            {'a', 'b', 'c', 'd',
                             'e', 'f', 'g', 'h',
                             'i', 'j', 'k', 'l',
                             'm', 'n', 'o', 'p',
                             'q', 'r', 's', 't',
                             'u', 'v', 'w', 'x',
                             'y', 'z'};

    unsigned int getCount(const string &str, unsigned int count, const vector<unsigned int> &suffixArray,
                          const vector<unsigned int> &lcp) const;
};


Solve::Solve() {
    readStr();
}


void Solve::readStr() {
    cin >> m_sString;

}


unsigned Solve::countDifferentSubstrings(string &str, vector<char> &alphabet) {
    unsigned count = 0;
    SuffArr sf(str, alphabet);
    vector<unsigned> suffixArray = sf.buildSuffixArray();
    vector<unsigned> lcp = sf.buildLCP();

    count = getCount(str, count, suffixArray, lcp);

    return count;
}


unsigned int Solve::getCount(const string &str, unsigned int count, const vector<unsigned int> &suffixArray,
                             const vector<unsigned int> &lcp) const {
    for (unsigned i = 0; i < str.size(); i ++)
        count += str.size() - suffixArray[i];

    for (unsigned i = 0; i < str.size() - 1; i ++)
        count -= lcp[i];
    return count;
}


void Solve::printAns() {
    cout << countDifferentSubstrings(m_sString, m_vcAlphabet) << endl;
}


//--------------------------------------------------------------------------------------------


int main () {
    Solve sl;
    sl.printAns();
    return 0;
}

//112 7,45
//110 7
//108 6,4
//106 6,35 6,45 106
//104 5,8 5,8 104
//102 5,7 5,7 102
//100 5,35 5,35 100
//98 4,95 4,95 98
//4,55 96
//
//95 4,35
//94 4,1 4,2 94
//92 3,8 3,8 92
//90 3,45 3,4 90
//88 3,1
//
//refactor 199-212(!)
//
//TODO:
//1)-
//2)--solve to alpha-
//3)/
//4)makeclass const(???)
