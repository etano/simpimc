#include "WriteClass.h"

void Writes::DoEvent()
{
  vector<Event*>::iterator iter;
  for (iter=events.begin(); iter!=events.end(); ++iter) {
    (*iter)->Write();
  }
}
