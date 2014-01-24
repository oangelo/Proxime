#ifndef TEST_FUNCTIONS_H
#define TEST_FUNCTIONS_H 
#include "../src/functions/rossler.h"

TEST(Functions, Rossler) {
    container variables({1, 1, 1});
    container parameters({1, 1, 1});
    value t;
    
    RosslerFunction func;
    func.set(t, variables, parameters);
    EXPECT_EQ(func.get_result(0), -2);
    EXPECT_EQ(func.get_result(1), 2);
    EXPECT_EQ(func.get_result(2), 1);
}

#endif /* TEST_FUNCTIONS_H */
