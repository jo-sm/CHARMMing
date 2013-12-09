import datetime
from qsar.models import attributes, object_attributes

class objQSARJob():
    id=0
    scheduler_id=0
    owner_index=0
    description=""
    type=""
    status=""
    queued=datetime.datetime.now()
    waited=datetime.datetime.now()
    runtime=datetime.datetime.now()
    totaltime=datetime.datetime.now()
    output_file=""
    output_data=""
    qsar_model_id=""
    target_conformations_file_list=[]
    ligands_file_list=[]
    results_file_list=[]

def AssignObjectAttribute(user_id, object_id, object_table, attribute_name, attribute_value):

    #logfile=open("/tmp/qsar_attr.log", "w")
    try:
        qsar_attribute=attributes.objects.get(attribute_short_name=attribute_name,attribute_owner=user_id)
        #logfile.write("attribute found %s\n" % attribute_name)
    except:
        qsar_attribute=attributes()
        qsar_attribute.attribute_short_name=attribute_name
        qsar_attribute.attribute_owner_id=user_id
        qsar_attribute.save()
        #logfile.write("new attribute created %s\n" % attribute_name)
    
    try:
    #if 1==1:
        qsar_object_attribute=object_attributes.objects.get( \
                              object_attribute_owner=user_id,attribute=qsar_attribute.id, \
                              object_table_name=object_table,object_id=object_id)
        qsar_object_attribute.attribute_value=attribute_value
        qsar_object_attribute.save()
        #logfile.write("updating object attribute %s\n")
    except:
    #else:
    #if 1==1:
        #logfile.write("creating new object attribute\n") 
        new_object_attribute=object_attributes()
        new_object_attribute.object_attribute_owner_id=user_id
        new_object_attribute.attribute_id=qsar_attribute.id
        new_object_attribute.object_table_name=object_table
        new_object_attribute.object_id=object_id
        new_object_attribute.value=attribute_value
        new_object_attribute.save()
    
    #logfile.close()

def GetObjectAttributeValue(user_id, object_id, object_table, attribute_name):
    #try:
    qsar_attribute=attributes.objects.get(attribute_short_name=attribute_name,attribute_owner=user_id)
    object_attribute=object_attributes.objects.get(object_attribute_owner=user_id,object_id=object_id, \
                     attribute=qsar_attribute,object_table_name=object_table)
    return object_attribute.value

def RemoveObjectAttributeValue(user_id, object_id, object_table, attribute_name):
    #try:                                                                                                 
    qsar_attribute=attributes.objects.get(attribute_short_name=attribute_name,attribute_owner_id=user_id)
    object_attribute=object_attributes.objects.get(object_attribute_owner=user_id,object_id=object_id, \
                     attribute=qsar_attribute,object_table_name=object_table)
    object_attribute.delete()
    #return object_attribute.value

class ObjectWithAttributes():
    
    object_id=""
    object_owner_id=""
    object_owner_index=""
    object_name=""
    object_type=""
    object_status=""
    attributes=[]

    def FillObjectAttributeList(self,object_table):
        self.attributes=object_attributes.objects.select_related().filter(object_id=self.object_id, object_table_name=object_table, object_attribute_owner=self.object_owner_id)
