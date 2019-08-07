extern "C" {
void ProfilerStart();

void gperfstart_() { ProfilerStart(); }
}
