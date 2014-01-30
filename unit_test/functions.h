#ifndef TEST_FUNCTIONS_H
#define TEST_FUNCTIONS_H 
#include "../src/functions/rossler.h"

TEST(Functions, get) {
    labels_values parameters{{"a", 1.},{"b", 1.},{"c", 1.},};
    container variables({1, 1, 1});
    value t;
    
    RosslerFunction func(parameters);
    func.set(t, variables);
    EXPECT_EQ(func["x"], -2);
    EXPECT_EQ(func["y"], 2);
    EXPECT_EQ(func["z"], 1);
}

TEST(Functions, Labels) {
    labels_values parameters{{"a", 1.},{"b", 1.},{"c", 1.},};
    container variables({1, 1, 1});
    value t;
    
    RosslerFunction func(parameters);
    func.set(t, variables);
    labels_values result(func.get_labels_values());

    EXPECT_EQ(result["x"], -2);
    EXPECT_EQ(result["y"], 2);
    EXPECT_EQ(result["z"], 1);
}


#endif /* TEST_FUNCTIONS_H */
