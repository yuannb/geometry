/*
 * sample http server
 * more detail see examples/httpd
 *
 */

#include "HttpServer.h"
#include "hthread.h"    // import hv_gettid
#include "hasync.h"     // import hv::async
#include <json.hpp>
#include <vector>
#ifdef blocking
    #undef blocking
#endif

#ifdef max 
    #undef max
#endif

#ifdef min 
    #undef min
#endif
#include "ReadJson.hpp"
#include "databaseManager.hpp"

using json = nlohmann::ordered_json;


using namespace hv;

/*
 * #define TEST_HTTPS 1
 *
 * @build   ./configure --with-openssl && make clean && make
 *
 * @server  bin/http_server_test 8080
 *
 * @client  curl -v http://127.0.0.1:8080/ping
 *          curl -v https://127.0.0.1:8443/ping --insecure
 *          bin/curl -v http://127.0.0.1:8080/ping
 *          bin/curl -v https://127.0.0.1:8443/ping
 *
 */
#define TEST_HTTPS 0

int main(int argc, char** argv) {

    int port = 0;
    if (argc > 1) {
        port = atoi(argv[1]);
    }
    if (port == 0) port = 8080;

    std::cout << port << std::endl;

    HttpService router;

    /* Static file service */
    // curl -v http://ip:port/
    router.Static("/", "./html");

    /* Forward proxy service */
    router.EnableForwardProxy();
    // curl -v http://httpbin.org/get --proxy http://127.0.0.1:8080
    router.AddTrustProxy("*httpbin.org");

    /* Reverse proxy service */
    // curl -v http://ip:port/httpbin/get
    router.Proxy("/httpbin/", "http://httpbin.org/");

    /* API handlers */
    // curl -v http://ip:port/ping
    router.GET("/ping", [](HttpRequest* req, HttpResponse* resp) {
        return resp->String("pong");
    });

    // curl -v http://ip:port/data
    router.GET("/data", [](HttpRequest* req, HttpResponse* resp) {
        static char data[] = "0123456789";
        return resp->Data(data, 10 /*, false */);
    });

    // curl -v http://ip:port/paths
    router.GET("/paths", [&router](HttpRequest* req, HttpResponse* resp) {
        return resp->Json(router.Paths());
    });

    // curl -v http://ip:port/get?env=1
    router.GET("/get", [](const HttpContextPtr& ctx) {
        hv::Json resp;
        resp["origin"] = ctx->ip();
        resp["url"] = ctx->url();
        resp["args"] = ctx->params();
        resp["headers"] = ctx->headers();
        return ctx->send(resp.dump(2));
    });

    // curl -v http://ip:port/echo -d "hello,world!"
    router.POST("/echo", [](const HttpContextPtr& ctx) {
        
        // NodeGraphDBManager ngdb2;
        // std::string gid("380eb874-aba2-43b4-b696-bbbdc601bd21");
        // bool flag = ngdb2.delete_nodegraph(gid);
        // 
        // return ctx->send("true", ctx->type());
        
        
        // NodeGraphDBManager ngdb2;
        // ngdb2.create_tables();
        // std::string gid("380eb874-aba2-43b4-b696-bbbdc601bd21");
        // NodeGraph ngt = ngdb2.get_nodegraph(gid);
        // for (auto it = ngt.m_curves.begin(); it != ngt.m_curves.end(); ++it)
        // {
        //     ngt.m_curves_disc.insert(std::make_pair(it->first, curve_disc(it->second)));
        // }
        // 
        // for (auto it = ngt.m_surfaces.begin(); it != ngt.m_surfaces.end(); ++it)
        // {
        //     ngt.m_surface_disc.insert(std::make_pair(it->first, surf_disc(it->second)));
        // }

        // json result = ngt.to_json();
        // 
        // std::string ret_str = result.dump();
        // std::cout << ret_str << std::endl;

        // return ctx->send(ret_str, ctx->type());
        std::string temp = ctx->body();
        // std::cout << temp << std::endl;
        json graph;
        try
        {
			graph = nlohmann::ordered_json::parse(temp);
        }
        catch (const nlohmann::json::parse_error& e)
        {
            std::cerr << "Json error: " << e.what() << std::endl;
        }
        // std::cout << "333333333333333333333333333333333333" << std::endl;
        // std::cout << graph << std::endl;
        // json graph = json::parse(temp);

        if (graph.contains("id") == false)
        {
			return ctx->send(std::string("graph is invalid"), ctx->type());
        }

        const json id = graph.at(ParamterType::id);
        if (id.is_string() == false)
        {
			return ctx->send(std::string("graph is invalid"), ctx->type());
        }
        NodeGraphDBManager ngdb;
        ngdb.create_tables();
        NodeGraph ng(id.get<std::string>());
        CreateNodeGraph(graph, ng);
        ng.calculate();
       
        // NodeGraphDBManager ngdb;
        // ngdb.create_tables();
		std::string ng_id = id.get<std::string>();
        std::string ng_date;
        if (ng.from_database == false)
        {
            ngdb.delete_nodegraph(ng.m_id);
			ngdb.add_nodegraph(ng);
            ng_date = ngdb.get_ng_date(ng_id);
        }
        else
        {
            ngdb.update_ng_time(ng_id);

            for (auto it = ng.changed_nodes.begin(); it != ng.changed_nodes.end(); ++it)
            {
                ngdb.delete_node(*it);
                ngdb.add_node(*ng.m_nodes[*it], ng_id);
            }
            
            for (auto it = ng.deleted_nodes.begin(); it != ng.deleted_nodes.end(); ++it)
            {
                ngdb.delete_node(*it);
            }
            ng_date = ngdb.get_ng_date(ng_id);
        }

        json result = ng.to_json();
        result["date"] = ng_date;
        
        std::string ret_str = result.dump();
        // std::cout << ret_str << std::endl;
        std::cout << "-------------------------end-----------------------" << std::endl;
        return ctx->send(ret_str, ctx->type());
    });

    // curl -v http://ip:port/user/123
    router.GET("/user/{id}", [](const HttpContextPtr& ctx) {
        hv::Json resp;
        resp["id"] = ctx->param("id");
        return ctx->send(resp.dump(2));
    });

    // curl -v http://ip:port/async
    router.GET("/async", [](const HttpRequestPtr& req, const HttpResponseWriterPtr& writer) {
        writer->Begin();
        writer->WriteHeader("X-Response-tid", hv_gettid());
        writer->WriteHeader("Content-Type", "text/plain");
        writer->WriteBody("This is an async response.\n");
        writer->End();
    });

    // middleware
    router.AllowCORS();
    router.Use([](HttpRequest* req, HttpResponse* resp) {
        resp->SetHeader("X-Request-tid", hv::to_string(hv_gettid()));
        return HTTP_STATUS_NEXT;
    });

    HttpServer server;
    server.service = &router;
    server.port = port;
#if TEST_HTTPS
    server.https_port = 8443;
    hssl_ctx_opt_t param;
    memset(&param, 0, sizeof(param));
    param.crt_file = "cert/server.crt";
    param.key_file = "cert/server.key";
    param.endpoint = HSSL_SERVER;
    if (server.newSslCtx(&param) != 0) {
        fprintf(stderr, "new SSL_CTX failed!\n");
        return -20;
    }
#endif

    // uncomment to test multi-processes
    // server.setProcessNum(4);
    // uncomment to test multi-threads
    // server.setThreadNum(4);

    server.start();

    // press Enter to stop
    while (getchar() != '\n');
    hv::async::cleanup();
    return 0;
}