#include <iostream>
#include <armadillo>

#include "b_ism.h"
#include "controlled-simulation.h"
#include "model-selection-strategy.h"
#include "mombf-bridge.h"
#include "exec_timer.h"

#define DEBUG_PRINTING // Print some additional debug info. Comment to disable.


// info
//
// Werden in jedem loop überschrieben
// s2_i
// b_i
// w_i
//
// Bleiben Konstant
//

int main() {
    arma::mat data = {
        {1, 1, -23.45266676, -2.30936055, -0.55970183, 1.34953562},
        {1, 2, -7.91095366, -1.46139144, 0.21277444, -0.47043324},
        {1, 3, -10.58970802, -0.06862945, -0.12904591, 1.80418798},
        {1, 4, -4.03629912, -0.23848014, 0.31444212, 0.18225964},
        {1, 5, 20.07450224, 0.54917267, 0.18824194, -2.76781820},
        {1, 6, 14.98383100, 1.04541866, 0.01999001, -1.02430597},
        {1, 7, 4.36945658, -0.37946788, 0.32571192, -1.17450487},
        {1, 8, -8.35794726, -0.22696702, -0.28827137, 1.35956658},
        {1, 9, 10.42538523, 0.09199853, -0.72593770, -1.00533029},
        {1, 10, -2.68719828, -0.85342818, -0.29560190, 1.31471504},
        {1, 11, 12.63427833, -0.50768911, -0.93493273, -0.60211238},
        {1, 12, 15.25046563, 1.12089459, -0.43084369, 0.55911127},
        {1, 13, 20.55473321, 0.58630229, -0.68998967, -0.40238799},
        {1, 14, 7.69004829, -1.27681905, 0.06944209, -0.88038903},
        {1, 15, 0.01908629, -0.70259846, 0.75558385, 0.55643298},
        {1, 16, 17.14391279, 1.69326930, 1.02409623, -0.06464430},
        {1, 17, 16.79348380, -0.27937351, -0.68313724, -1.14558556},
        {1, 18, 18.42294041, 1.23860288, 0.89303727, -1.17465964},
        {1, 19, 15.93351524, -0.27851469, -1.51865248, -0.30728423},
        {1, 20, 3.41723235, 0.50750409, 2.13166805, 0.07164625},
        {2, 1, 2.01367442, -0.52964157, -0.84032731, -0.24943370},
        {2, 2, 0.12292821, 1.24425732, 1.89566249, -0.07214188},
        {2, 3, 14.91203433, 0.74228208, 0.57532928, -2.24319079},
        {2, 4, -5.11301574, -0.46578083, 0.89250538, -0.64141230},
        {2, 5, -2.93105027, 0.75886611, 0.61983879, 0.82792110},
        {2, 6, -1.85769585, 0.45025558, 0.09238917, 0.81542852},
        {2, 7, 12.58391461, 0.61893994, -0.07080554, -1.25784309},
        {2, 8, 6.08416188, 0.43634622, -0.35005450, -0.11806201},
        {2, 9, 4.04203852, -0.06634003, -0.35659074, -0.11128336},
        {2, 10, -5.61360271, -0.24980976, -0.11055986, 0.75599906},
        {2, 11, 6.40934135, 1.33751293, -0.26408725, 0.72592710},
        {2, 12, -11.60791505, -1.09736930, -0.61361462, 1.08042091},
        {2, 13, 4.77833559, 1.24434837, 1.19657744, -0.14257870},
        {2, 14, 3.98187096, 0.82541432, -0.90584227, 1.00333061},
        {2, 15, 4.86749764, 1.11334388, 0.66078066, -0.18836024},
        {2, 16, -2.97304347, 0.17946976, 0.54460379, 0.19905417},
        {2, 17, -1.24581570, -1.49029084, -0.48200822, -1.30089775},
        {2, 18, 10.14556265, 0.93577627, -0.38133499, -0.41726513},
        {2, 19, 8.77631832, 0.17082548, -0.53630861, -0.60742406},
        {2, 20, -5.32484498, -0.27814329, -0.17968005, 0.70029078},
        {3, 1, 6.34106622, 1.03389672, -1.18754295, 1.44397345},
        {3, 2, -0.10247840, -0.52907998, 0.71449362, -1.34239109},
        {3, 3, 0.93901800, 0.30634010, -0.78423320, 0.78976343},
        {3, 4, -8.23697372, -1.01782854, 0.22088162, 0.11148258},
        {3, 5, -9.49952944, -0.33937959, 0.82135615, 0.33237310},
        {3, 6, 3.24726127, -0.15963926, -0.69482875, -0.23291455},
        {3, 7, -8.94658045, -2.23058885, -0.09104611, -1.00355918},
        {3, 8, -18.67803117, -0.02682357, 1.79855827, 1.53172391},
        {3, 9, 18.56838167, 2.10492018, 1.05464263, -1.08119123},
        {3, 10, -12.34777296, -0.23481469, 2.47325089, 0.07734140},
        {3, 11, -19.62430162, -1.70293765, 1.14583713, 0.49988338},
        {3, 12, 4.42893308, 1.42905860, -0.12474232, 0.96249562},
        {3, 13, -4.50378125, -0.77566028, -0.58205656, 0.42390948},
        {3, 14, 7.53966294, 0.81322432, 0.35094792, -0.50595807},
        {3, 15, 14.47285868, 1.32722746, -2.44584316, 1.26806703},
        {3, 16, -1.47646208, -0.35883300, -0.44026205, 0.07695551},
        {3, 17, 6.69367344, 1.42208596, 1.30404222, -0.83074969},
        {3, 18, 9.61231686, 0.06315024, -1.24184881, -0.52511487},
        {3, 19, 10.14149005, 0.44613682, -0.63743304, -0.55301099},
        {3, 20, -0.69086139, -2.07200545, -1.89931838, -0.63308381}
    };

    // Assuming 'data' is your data matrix
    BIsmOutput result = b_ism(
        data,
        0,
        1,
        2,
        5000,
        500,
        "g",
        100.0,
        0.001,
        0.001,
        1.0,
        1.0,
        1.0,
        false,
        true,
        false,
        false,
        false,
        true,
        true
    );

    // Print the results
    std::cout << "Column means of w_store:" << std::endl;
    std::cout << result.w_store_means.t() << std::endl;


    // Print the results
    std::cout << "Column means of b_store:" << std::endl;
    std::cout << result.b_store_means.t() << std::endl;

    // Print the results
    std::cout << "Column means of s2_store:" << std::endl;
    std::cout << result.s2_store_means.t() << std::endl;

    return 0;
}


