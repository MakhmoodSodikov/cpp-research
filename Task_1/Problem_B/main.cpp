#include <vector>
#include <algorithm>
#include <string>
#include <iostream>
#include <fstream>


class Transformations {
public:
    std::vector<int> stringToZFunc(std::string& Text);
    std::string zFuncToString(std::vector<int>& ZFunc);
    std::vector<int> stringToPrefixFunc(std::string& Text);
    std::string prefixFuncToString(std::vector<int>& PrefixFunc);
    std::vector<int> prefixFuncToZFunc(std::vector<int>& PrefixFunc);
    std::vector<int> zFuncToPrefixFunc(std::vector<int> &ZFunc);
private:
    char generateSample(std::vector<int>& Pref, int i, std::string& text);
    std::vector<int> zFuncToPrefix(std::vector<int> &ZFunc);
};


std::vector<int> Transformations::stringToZFunc(std::string& Text) {
    std::vector<int> z_func(Text.length(), 0);
    z_func[0] = Text.length();

    int left_iter = 0, right_iter = 0;

    for (int iter = 1; iter < Text.length(); ++iter) {

        if (iter <= right_iter) {
            z_func[iter] = std::min(right_iter - iter + 1, z_func[iter - left_iter]);

            while (iter + z_func[iter] < Text.length() && Text[z_func[iter]] == Text[iter + z_func[iter]]) {
                z_func[iter]++;
            }
        } else {
            while (iter + z_func[iter] < Text.length() && Text[z_func[iter]] == Text[iter + z_func[iter]]) {
                z_func[iter]++;
            }
        }

        if (iter + z_func[iter] - 1 > right_iter) {
            right_iter = iter + z_func[iter] - 1;
            left_iter = iter;
        }
    }

    return z_func;
}


char Transformations::generateSample(std::vector<int>& Pref, int i, std::string& text) {
    int curr_pref = Pref[i - 1];
    char curr_sample = 'b';

    while (curr_pref != 0) {

        curr_sample = std::max(curr_sample, (char) (text[curr_pref] + 1));
        curr_pref = Pref[curr_pref - 1];

    }

    return curr_sample;
}


std::string Transformations::prefixFuncToString(std::vector<int>& PrefixFunc) {
    std::string text = "a";

    for (int i = 1; i < PrefixFunc.size(); ++i) {

        if (PrefixFunc[i] != 0) {
            text += text[PrefixFunc[i] - 1];
        } else {
            text += generateSample(PrefixFunc, i, text);
        }

    }

    return text;
}


std::vector<int> Transformations::zFuncToPrefix(std::vector<int> &ZFunc) {
    std::vector<int> prefix_func(ZFunc.size(), 0);

    for (int i = 1; i < ZFunc.size(); i ++) {
        if (ZFunc[i] != 0) {
            for (int j = ZFunc[i]; j >= 0 && prefix_func[i + j - 1] == 0; j --) {
                prefix_func[i + j - 1] = j;
            }
        }
    }

    return prefix_func;
}


std::string Transformations::zFuncToString(std::vector<int> &ZFunc) {
    std::vector<int> pref_func = zFuncToPrefix(ZFunc);
    std::string text = prefixFuncToString(pref_func);

    return text;
}


std::vector<int> Transformations::stringToPrefixFunc(std::string &Text) {
    std::vector<int> prefix_func(Text.length(), 0);

    for (int iter = 1; iter < Text.length(); ++iter) {
        int sample_iter = prefix_func[iter - 1];

        if (Text[iter] == Text[sample_iter]) {
            prefix_func[iter] = sample_iter + 1;
        } else {
            while (Text[sample_iter] != Text[iter] && sample_iter > 0) {
                sample_iter = prefix_func[sample_iter - 1];
            }

            if (Text[iter] == Text[sample_iter]) {
                prefix_func[iter] = sample_iter + 1;
            } else {
                if (sample_iter == 0) {
                    prefix_func[iter] = 0;
                }
            }
        }
    }

    return prefix_func;
}


std::vector<int> Transformations::prefixFuncToZFunc(std::vector<int> &PrefixFunc) {
    std::string text = prefixFuncToString(PrefixFunc);

    std::vector<int> z_func = stringToZFunc(text);

    return z_func;
}


std::vector<int> Transformations::zFuncToPrefixFunc(std::vector<int> &ZFunc) {
    std::string text = zFuncToString(ZFunc);

    std::vector<int> pref_func = stringToPrefixFunc(text);

    return pref_func;
}


void problem_B1(Transformations transformations) {
    std::vector<int> pref;

    std::ifstream file("input.txt");

    while (!file.eof()) {
        int num;
        file>>num;
        pref.push_back(num);
    }

    std::cout << transformations.prefixFuncToString(pref);

    file.close();
}


void problem_B2(Transformations transformations) {
    std::vector<int> z_func;
    int n;

    while (std::cin >> n)
        z_func.push_back(n);

    std::cout << transformations.zFuncToString(z_func);
}


std::ostream& operator << (std::ostream& os, const std::vector<int>& vec) {
    for (auto elem : vec) {
        std::cout << elem << ' ';
    }
    return os;
}


void tester(Transformations transformations) {
    std::string str = "abacaba";
    std::vector<int> z_func = transformations.stringToPrefixFunc(str);
    std::vector<int> pref_func = transformations.stringToZFunc(str);
    //----------
    std::cout << str;
    std::cout << z_func;
    std::cout << pref_func;
    /*
     * tests:
                        'a b a c a b a'
    'abacaba' to z     : 7 0 1 0 3 0 1
    'abacaba' to prefix: 0 0 1 0 1 2 3
    */
}


int main() {
    Transformations tr;
//    problem_B1(tr);
//    problem_B2(tr);
//    tester(tr);
    return 0;
}
