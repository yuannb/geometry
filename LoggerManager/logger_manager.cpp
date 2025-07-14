#include "logger_manager.hpp"
#include "spdlog/sinks/basic_file_sink.h"
#include "spdlog/sinks/stdout_color_sinks.h"

using namespace spdlog::sinks;
std::shared_ptr<spdlog::logger> LoggerManager::getLogger()
{
    static std::shared_ptr<spdlog::logger> logger;
    if (!logger)
    {
		const std::string log_name = "multi_sink";
        // 创建文件 sink，设置日志文件路径
        auto file_sink = std::make_shared<basic_file_sink_mt>(log_name, true);
		file_sink->set_level(spdlog::level::info);

		auto console_sink = std::make_shared<stdout_color_sink_mt>();
		console_sink->set_level(spdlog::level::err);

		spdlog::sinks_init_list sink_list = { console_sink, file_sink };
		logger = std::make_shared<spdlog::logger>(log_name, sink_list);
		logger->set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%n] [%l] %v");
		logger->set_level(spdlog::level::info);

        // 注册到 spdlog 全局注册表，方便库内部获取
        spdlog::register_logger(logger);
    }
    return logger;
}