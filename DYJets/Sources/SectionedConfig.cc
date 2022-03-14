#include "SectionedConfig.h"

template <> int SectionedConfig::get(const char *section, const char *param) const
{
    return get(section, param, 0);
}

template <> unsigned SectionedConfig::get(const char *section, const char *param) const
{
    return get(section, param, 0);
}

template <> bool SectionedConfig::get(const char *section, const char *param) const
{
    return get(section, param, false);
}

template <> float SectionedConfig::get(const char *section, const char *param) const
{
    return get(section, param, false);
}

template <> double SectionedConfig::get(const char *section, const char *param) const
{
    return get(section, param, false);
}
