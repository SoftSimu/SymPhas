
#include "testnextsave.h"


void testnextsave()
{
	iter_type saveme = symphas::conf::config().save.get_start();
	auto save = symphas::conf::config().save;

	// this test assumes that the base values give the correct next value
	for (iter_type i = symphas::conf::config().save.get_start(); i <= symphas::conf::config().save.get_stop(); i = symphas::conf::next_save(i), ++saveme)
	{
		iter_type next = symphas::conf::next_save(i);
		printf("\ntesting from %d to %d\n...", i, next);
		for (iter_type j = i; j < next; ++j)
		{
			if (symphas::conf::next_save(j) != next)
			{
				printf("failed on index = %d; should have equalled %d but was %d\n", 
					j, 
					next, 
					symphas::conf::next_save(j));
			}
		}
	}

	if (saveme - 1 != symphas::conf::config().save.num_saves())
	{
		printf("the number of saves are not consistent with what's computed\n");
	}

}

