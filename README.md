# 照护衔接与等待分析

本项目使用 Python 3.11 及以上版本，仅依赖标准库，无需安装第三方包。数据文件均为本地合成记录，不含真实个人或机构资料。

encounters 校验机构、人员、住院和转介，保存修订；pathways 构建明确关联的照护链；waiting_analysis 分析等待曲线及计划；exports 导出内部和公开结果。

仅使用明确的身份及转介关联。业务截止与知识截止独立保存。现有等待曲线排除拒绝和取消；现有住院过程曲线将任何出院计为事件。统计描述观测记录，不提供临床建议。

在项目根目录运行完整测试：

```sh
python3 -m unittest discover -s tests -v
```

编译检查：

```sh
python3 -m compileall -q .
```

实际演示：

```sh
python3 continuity.py demo --fixture fixtures/care.json
```

演示在临时目录内建立本地状态，结束后自动清理。测试、编译与演示在依赖已准备的环境内均不下载数据或软件。输入错误以非零退出码及错误对象报告。更多子命令可通过对应入口的 --help 查看。

容器环境可执行 `docker build -t local-analysis .`，再运行 `docker run --rm --network none local-analysis`。镜像固定基础内容摘要，构建后运行只访问本地文件。完整测试与演示也可直接使用已有 Python 环境。
