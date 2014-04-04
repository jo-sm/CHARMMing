--
-- Dumping data for table `qsar_job_types`
--

LOCK TABLES `qsar_job_types` WRITE;
/*!40000 ALTER TABLE `qsar_job_types` DISABLE KEYS */;
INSERT INTO `qsar_job_types` VALUES (1,'Train',NULL),(2,'Predict',NULL);
/*!40000 ALTER TABLE `qsar_job_types` ENABLE KEYS */;
UNLOCK TABLES;

--
-- Dumping data for table `qsar_model_types`
--

LOCK TABLES `qsar_model_types` WRITE;
/*!40000 ALTER TABLE `qsar_model_types` DISABLE KEYS */;
INSERT INTO `qsar_model_types` VALUES (1,'Random Forest (QSAR) Regression',''),(2,'Random Forest (SAR) Categorization',''),(3,'SVM (SAR) Categorization',NULL),(4,'SVM (QSAR) Regression',NULL),(5,'Logit (SAR) Categorization',NULL),(6,'SGD (SAR) Categorization',NULL),(7,'Naive Bayes (SAR) Categorization',NULL),(8,'Gradient Boosting (SAR) Categorization',NULL),(10,'Decision Tree (SAR) Categorization',NULL),(11,'Gradient Boosting (QSAR) Regression',NULL),(12,'Decision Tree (QSAR) Regression',NULL),(14,'Ridge (QSAR) Regression',NULL),(15,'Lasso (QSAR) Regression',NULL),(16,'Elastic Net (QSAR) Regression',NULL),(17,'SGD (QSAR) Regression',NULL);
/*!40000 ALTER TABLE `qsar_model_types` ENABLE KEYS */;
UNLOCK TABLES;
