#pragma once
#include <spdlog/spdlog.h>
#include <memory>

class LoggerManager {
public:
    static std::shared_ptr<spdlog::logger> getLogger();
};