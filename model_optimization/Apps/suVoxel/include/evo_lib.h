#include <vector>
#include <iostream>
class rule_function
{
public:
	virtual void evolution()
	{
		std::cout << "This is the father class of the evolution rule";
	}
};	
class rule1 :public rule_function
{
public:
	void evolution()
	{
		cout << "sample";
	}
	static rule_function* gonly_;
	static rule_function* gonly();

};
rule_function* rule1::gonly_ = 0; 
rule_function* rule1::gonly()
{
	if (!gonly_)
	{
		gonly_ = new rule1;
		return gonly_;
	}
}


class rule_manager
{
public:
	std::vector<rule_function *> rule_list;
	void add_rule(rule_function *);

};
void rule_manager::add_rule(rule_function *r)
{
	rule_list.push_back(r);
}