from skavengers import Pipeline, Executor

if __name__ == '__main__':
    process_executor = Executor(3)
    pipeline_config = './skavengers/pipeline_setup.ini'
    pipeline = Pipeline(None, pipeline_config)
    pipeline.init()
    pipeline.execute(process_executor)
